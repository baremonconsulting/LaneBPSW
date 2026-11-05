# LaneBPSW — nativní validace a benchmark na Azure

Datum: 2026-11-04 (UTC)
Výsledek: **FUNCTIONAL PASS / BENCHMARK PASS**

## Rozhodnutí

LaneBPSW byl v tomto doloženém běhu sestaven a spuštěn nativně bez Intel SDE na procesoru, který přímo poskytuje AVX-512F a AVX-512 IFMA. V provedeném rozsahu nebyla zjištěna funkční neshoda. Výkonový běh ukazuje přibližně 4,0× až 4,8× zrychlení zahřáté operace `nextprime` podle průměrného thread CPU času proti systémovému GMP 6.3.0 a podle velikosti vstupu.

Výsledek potvrzuje použitelnost implementace na této konkrétní Azure/AMD konfiguraci. Není univerzální garancí výkonu na jiných procesorech nebo VM.

## Prostředí

- Azure VM, host `PGTESTING`, 4 vCPU;
- AMD EPYC 9V74 80-Core Processor, Microsoft full virtualization;
- dostupné příznaky: `avx512f`, `avx512ifma`, `avx512bw`, `avx512dq`, `avx512vl`;
- Ubuntu 24.04.5 LTS;
- Clang 18.1.3, `-O3 -march=native`, C++17;
- systémové GMP 6.3.0 (`2:6.3.0+dfsg-2ubuntu6.1`);
- OpenSSL 3.0.13 vývojové knihovny;
- hlavní benchmark připnutý na vCPU 1;
- zdrojový archiv LaneBPSW: SHA-256 `dcefb27fd8d14935d35b184ae758a0aff01191be448cfda7ea984c0e9dfac956`.

Úplný výpis `lscpu`, verzí, kernelu a hashů zdrojů je v `environment.txt`.

## Funkční validace

Nativní `make test` skončil návratovým kódem 0.

| Test | Kontroly | Neshody | Stav |
|---|---:|---:|---|
| `test_kernels` | 187 197 | 0 | PASS |
| `test_api` | 73 959 | 0 | PASS |
| `test_lucv` | 8 | 0 | PASS |
| `test_paths` | 1 926 | 0 | PASS |
| `test_bpsw --checks-only` | 5 201 | 0 | PASS |
| **Celkem** | **268 291** | **0** | **PASS** |

Dodatečný benchmarkový běh `test_bpsw` znovu provedl 5 201 kontrol bez neshody. Test chiavdf predikátu provedl dalších 2 105 porovnání lane-packed a scalar varianty: 0 neshod.

## Celá operace nextprime

Každý dotaz používá shodný vstup pro tři varianty a rotuje jejich pořadí:

1. GMP `mpz_nextprime`;
2. lane filtr + GMP finální test;
3. lane filtr + finální lane-packed BPSW.

Všechny varianty v každém řádku vrátily shodný výsledek.

| Bitů | Dotazů | GMP mean | Finální LaneBPSW mean | Zrychlení vs GMP | Zrychlení vs starší lane+GMP final |
|---:|---:|---:|---:|---:|---:|
| 256 | 3 000 | 0,1710 ms | 0,0432 ms | 3,96× | 2,00× |
| 512 | 1 500 | 1,2551 ms | 0,2606 ms | 4,82× | 1,77× |
| 768 | 800 | 4,1493 ms | 0,9603 ms | 4,32× | 1,50× |
| 1 024 | 600 | 11,5545 ms | 2,6001 ms | 4,44× | 1,39× |
| 1 536 | 200 | 50,4699 ms | 11,7416 ms | 4,30× | 1,25× |
| 2 048 | 150 | 138,0625 ms | 30,5315 ms | 4,52× | 1,22× |
| 3 072 | 60 | 619,9677 ms | 133,2452 ms | 4,65× | 1,15× |
| 4 096 | 40 | 1 780,9011 ms | 372,5228 ms | 4,78× | 1,11× |

Součet mismatchů: **0**. BPSW fallbacky v této benchmarkové sadě: **0**.

## Samotný finální primality test

Medián lane-packed BPSW proti `mpz_millerrabin(n,25)` na prvočíslech:

| Bitů | Lane median | GMP median | Poměr GMP/Lane |
|---:|---:|---:|---:|
| 256 | 17,1 µs | 59,7 µs | 3,49× |
| 512 | 76,1 µs | 259,3 µs | 3,41× |
| 768 | 266,3 µs | 697,8 µs | 2,62× |
| 1 024 | 530,4 µs | 1 601,2 µs | 3,02× |
| 1 536 | 1 777,3 µs | 4 891,7 µs | 2,75× |
| 2 048 | 4 049,5 µs | 11 100,2 µs | 2,74× |
| 3 072 | 13 082,2 µs | 33 414,9 µs | 2,55× |
| 4 096 | 29 720,5 µs | 69 236,8 µs | 2,33× |

## Interpretace

- Emulační funkční výsledek byl potvrzen na skutečných AVX-512 IFMA instrukcích.
- Kód funguje i na AMD EPYC, nikoliv pouze na dříve emulovaném či dokumentovaném Intel procesoru.
- Největší přínos celé `nextprime` operace proti GMP je v tomto běhu stabilně přibližně čtyř- až téměř pětinásobný.
- Přínos lane-packed finálního BPSW proti starší variantě s GMP final testem klesá s velikostí z přibližně 2,00× na 1,11×, protože u velkých vstupů větší část celkového času připadá na jiné fáze hledání.
- Absolutní časy jsou výsledkem jedné konkrétní virtualizované Azure instance. Pro tvrzení o jiném hardwaru je nutné měření zopakovat.

## Omezení

- VM je virtualizovaná a může být ovlivněna hostitelským plánováním a kolísáním frekvence.
- Hlavní běh obsahuje mnoho párových dotazů, ale není to více nezávislých běhů v různých časech/instancích.
- Nešlo o dlouhodobý vícevláknový test; knihovna sama není thread-safe.
- Implementace není constant-time a je určena pro veřejné vstupy.
- Výsledek neověřuje jiné verze GMP než 6.3.0 ani jiné ABI/architektury.

## Důkazní soubory

- `functional-tests.log`, `functional-tests.stderr`, `functional-tests.status`;
- `benchmark/np3.log` a osm CSV se všemi jednotlivými časy;
- `benchmark/test_bpsw.log`;
- `benchmark/bench_pass.log`;
- `benchmark/prof_np2.log`;
- `benchmark/hp_check.log`;
- `environment.txt`, `benchmark/env.txt`;
- `run_all.console.log`, `run_all.stderr`, `run_all.status`.
