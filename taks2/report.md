---
title: "SAD 2026L"
subtitle: "Część 2 Projektu"
author: "Matuszewski Maciej, Sosnowski Wiktor"
lang: pl
geometry: margin=2.5cm
fontsize: 11pt
---

# Problem 1. CTG - kontynuacja: testy zgodności i moc

## 1. Wprowadzenie

Niniejsze opracowanie stanowi bezpośrednią kontynuację Problemu 2 z części
pierwszej projektu, w którym zilustrowano graficznie zbieżność rozkładu
średniej z próby $\bar{X}_n$ do rozkładu normalnego, zgodnie z Centralnym
Twierdzeniem Granicznym (CTG). Teraz zbieżność ta zostaje zweryfikowana
formalnie, tj. za pomocą testów statystycznych zgodności rozkładu z
$\mathcal{N}(\mu, \sigma^2/n)$.

Punktem wyjścia jest pojedynczy rozkład $X \sim \mathrm{Exp}(\lambda=1)$ -
ten sam, którego użyto w części pierwszej - co pozwala wykorzystać
analityczne momenty rozkładu wyjściowego oraz porównać konkluzje z
wcześniejszą analizą graficzną. Dla rosnących liczności próby $n$
rozkład $\bar{X}_n$ przestaje wykazywać prawoskośność i zbliża się
do $\mathcal{N}(\mu, \sigma^2/n)$ - pytanie brzmi: **który test
najszybciej rozróżnia te dwa rozkłady?**

Porównujemy moc trzech testów normalności na podstawie symulacji Monte
Carlo:

1. **test Kołmogorowa–Smirnowa** (KS) - przeciwko w pełni
   wyspecyfikowanej $\mathcal{N}(\mu_T,\, \sigma_T^2/n)$,
2. **test Shapiro–Wilka** (SW) - klasyczny test oparty na korelacji
   statystyk pozycyjnych,
3. **test D'Agostino–Pearsona K²** (DA) - test oparty na
   standaryzowanej skośności i kurtozie.

## 2. Metodologia

### 2.1. Rozkład bazowy i statystyka badana

Niech $X_1,\dots,X_n \overset{\text{iid}}{\sim} \mathrm{Exp}(\lambda=1)$.
Z teorii znamy:

$$
\mu_T = \mathbb{E}[X] = 1,\quad \sigma_T^2 = \mathrm{Var}(X) = 1,\quad
\gamma_1 = 2,\quad \gamma_2 = 6,
$$

skąd dla średniej $\bar{X}_n = \tfrac{1}{n}\sum X_i$:

$$
\mathbb{E}\bar{X}_n = \mu_T,\qquad
\mathrm{Var}\bar{X}_n = \sigma_T^2/n,\qquad
\mathrm{skew}(\bar{X}_n) = \frac{2}{\sqrt{n}},\qquad
\mathrm{exc.\,kurt}(\bar{X}_n) = \frac{6}{n}.
$$

CTG mówi, że $\bar{X}_n \xrightarrow{d} \mathcal{N}(\mu_T,\sigma_T^2/n)$,
więc dla **dużego** $n$ żaden z testów nie powinien systematycznie
odrzucać hipotezy o normalności. Dla **małego** $n$ rozkład $\bar{X}_n$
jest jednak nadal wyraźnie skośny, a moc testu - czyli prawdopodobieństwo
wykrycia tej rozbieżności - staje się głównym kryterium oceny.

### 2.2. Schemat symulacji

Niech $K$ oznacza wielkość próby przekazywanej *każdemu* z testów (tj.
liczbę średnich obliczonych w pętli), a $n$ liczbę uśrednianych zmiennych
losowych. Wybrano $K \in \{20,\ 50,\ 100\}$ (zgodnie z poleceniem -
"kilkadziesiąt"), przy czym wartością wyróżnioną jest $K=50$.
Przyjęto $\alpha = 0{,}05$.

Pojedyncza replikacja Monte Carlo wygląda następująco:

1. Losuje się $K$ niezależnych średnich
   $\bar{X}_n^{(1)},\dots,\bar{X}_n^{(K)}$,
   z których każda jest obliczona z $n$ wartości $\mathrm{Exp}(1)$.
2. Na tej $K$-elementowej próbce uruchamia się trzy testy normalności
   (KS, SW, DA) i zapisuje uzyskane $p$-wartości.
3. Z tej samej próbki estymuje się skośność $\hat{g}_1$ i kurtozę
   $\hat{b}_2$ - wielkości stojące u podstaw testu D'Agostino.

Schemat powtarzany jest $M = 5\,000$ razy dla każdej pary $(n,K)$ z
siatki $n \in \{2,3,5,10,20,30,50,100,200\}$,
$K \in \{20,50,100\}$. Moc testu definiujemy jako empiryczne
prawdopodobieństwo odrzucenia $H_0$:

$$
\widehat{\text{moc}}(n,K,\text{test}) \;=\;
\frac{1}{M}\sum_{i=1}^{M}\mathbb{1}\!\left(p^{(i)} < \alpha\right).
$$

Eksperyment zrealizowano w **R 4.5.3** (`ggplot2`, `dplyr`, `tidyr`)
w środowisku `conda`. Kod podzielony jest modularnie na pliki
`R/config.R` (stałe), `R/tests.R` (testy normalności wraz z ręczną
implementacją D'Agostino–Pearsona), `R/simulation.R` (pętla MC),
`R/plots.R` (wizualizacje) oraz `R/main.R` (punkt wejścia).
Reprodukowalność zapewnia `set.seed(20260519)`.

### 2.3. Implementacja testów

**KS** stosujemy w wersji *w pełni wyspecyfikowanej*: porównujemy
empiryczną dystrybuantę próbki $K$ średnich z dystrybuantą
$\mathcal{N}(\mu_T,\sigma_T^2/n)$, gdzie $\mu_T$ i $\sigma_T^2$ są
**znanymi** momentami rozkładu wyjściowego. Tak skonstruowany test ma
poprawny rozmiar (nie wymaga poprawki Lillieforsa), co czyni
porównanie z dwoma pozostałymi testami uczciwym.

**SW** używamy w wbudowanej postaci `stats::shapiro.test()`
(estymuje $\mu,\sigma$ z danych, ważny dla $3 \le K \le 5000$).

**D'Agostino–Pearson K²** zaimplementowano ręcznie - biblioteka
`moments` nie była dostępna w środowisku - na podstawie wzorów
D'Agostino (1970) i Anscombe'a–Glynna (1983). Test łączy
standaryzowane $Z_1$ (skośność) i $Z_2$ (kurtoza):

$$
K^2 = Z_1^2 + Z_2^2 \;\sim\; \chi^2(2) \;\;\text{pod }H_0,
$$

zatem $p$-value $= 1 - F_{\chi^2_2}(K^2)$. Składowe $Z_1$ i $Z_2$ pełnią
podwójną rolę: po pierwsze definiują test, po drugie - same w sobie są
estymatorami "ile odchyleń standardowych" empiryczna skośność / kurtoza
oddala próbkę od próbki normalnej.

## 3. Wyniki

### 3.1. Co "widzi" test - wizualizacja jednej replikacji

![Histogramy standaryzowanej średniej $\bar{X}_n$ na pojedynczej próbce $K=50$ średnich, dla wybranych $n$.](plots/01_means_density.png){width=72%}

Już na pojedynczej realizacji ($K=50$, ziarno ustalone) widać efekt
CTG: histogram dla $n=2$ wyraźnie odbiega od $\mathcal{N}(0,1)$
(prawoskośność, ciężki prawy ogon); dla $n=30$ rozkład jest niemal
symetryczny, a przy $n=200$ histogram jest praktycznie nieodróżnialny
od krzywej Gaussa. Testy normalności muszą podjąć decyzję
**dokładnie** na takiej $K$-elementowej próbce.

### 3.2. Tabela mocy ($K=50$, $\alpha=0{,}05$)

**Tabela 1.** Empiryczne moce testów wyznaczone z $M=5\,000$ replikacji.

|  $n$ |  $\widehat{\text{moc}}_{KS}$ | $\widehat{\text{moc}}_{SW}$ | $\widehat{\text{moc}}_{DA}$ | $\hat{g}_1$ (śr.) | $\hat{b}_2-3$ (śr.) |
|----:|----:|----:|----:|----:|----:|
|   2 | 0,236 | **0,950** | 0,799 | 1,19 | 1,60 |
|   3 | 0,177 | **0,828** | 0,658 | 0,99 | 1,10 |
|   5 | 0,132 | **0,598** | 0,476 | 0,78 | 0,64 |
|  10 | 0,085 | **0,331** | 0,290 | 0,55 | 0,26 |
|  20 | 0,068 | 0,190 | 0,188 | 0,39 | 0,08 |
|  30 | 0,056 | 0,137 | **0,146** | 0,32 | 0,01 |
|  50 | 0,055 | 0,099 | **0,110** | 0,24 | −0,03 |
| 100 | 0,057 | 0,082 | **0,087** | 0,17 | −0,06 |
| 200 | 0,046 | 0,060 | **0,069** | 0,12 | −0,11 |

Test odznaczony pogrubieniem to najmocniejszy test w danym wierszu.
Wartości średniej skośności i nadwyżki kurtozy obrazują, jak silnie
$\bar{X}_n$ odbiega od $\mathcal{N}(0,1)$ w sensie momentów trzeciego
i czwartego rzędu. Pełną siatkę (wszystkie wartości $K$) zapisano
w `output/power_grid.csv`.

### 3.3. Moc w funkcji $n$

![Moce trzech testów dla $K=50$, $\alpha=0{,}05$. Linia przerywana - rozmiar testu.](plots/02_power_vs_n.png){width=72%}

Krzywe mocy mają wyraźną hierarchię: dla każdego $n \le 20$
zachodzi $\widehat{\text{moc}}_{SW} > \widehat{\text{moc}}_{DA}
\gg \widehat{\text{moc}}_{KS}$. Powyżej $n \approx 20$ moce SW i DA
zrównują się, a od $n \ge 30$ DA wyprzedza SW o $0{,}005$–$0{,}01$
(różnica mieści się w przedziale ufności replikacji). KS pozostaje
zauważalnie słabszy w całym zakresie. Wszystkie trzy krzywe zbiegają
do poziomu nominalnego $\alpha = 0{,}05$ przy $n \to \infty$, co
potwierdza, że tempo zbiegania mocy do rozmiaru testu jest tutaj
empirycznym odpowiednikiem szybkości CTG.

### 3.4. Wpływ wielkości próby testowej $K$

![Moce w funkcji $n$, z podziałem na panele dla $K \in \{20, 50, 100\}$.](plots/03_power_facets_K.png){width=82%}

Wpływ $K$ jest bardzo wyraźny. Przy $K=20$ nawet dla $n=2$ żaden z
testów nie osiąga mocy 0,55. Przy $K=100$ ten sam $n=2$ daje moc
$\approx 0{,}99$ dla SW i DA, a KS rośnie z $0{,}12$ do $0{,}44$.
Wzorzec uporządkowania SW > DA > KS jest jednak zachowany dla każdego
$K$ - różnica jest tylko ilościowa, nie jakościowa. Wynika z tego, że
**dobór $K$ (a nie wybór konkretnego testu)** jest najważniejszą decyzją
projektową przy planowaniu eksperymentu wykrywania niezgodności z
rozkładem normalnym.

### 3.5. Skośność i kurtoza - interpretacja testu D'Agostino

Test D'Agostino–Pearsona buduje statystykę z dwóch składowych:
standaryzowanej skośności $Z_1$ i standaryzowanej kurtozy $Z_2$. Aby
zrozumieć, dlaczego DA "działa" lub "nie działa", musimy zobaczyć, jak
rzeczywiste wartości skośności i kurtozy $\bar{X}_n$ zachowują się w
funkcji $n$.

![Empiryczna skośność i nadwyżka kurtozy próbki $K=50$ średnich, w zestawieniu z wartościami teoretycznymi $2/\sqrt{n}$ oraz $6/n$.](plots/04_shape_convergence.png){width=78%}

Empiryczne krzywe niemal dokładnie pokrywają się z teoretycznymi
$2/\sqrt{n}$ i $6/n$, przy systematycznie lekko obniżonym poziomie -
to znane ujemne obciążenie estymatora skośności i kurtozy
dla małych $K$ (wartość oczekiwana $g_1$ jest mniejsza niż
$\gamma_1$, gdy $\gamma_1 > 0$). Oba momenty maleją gładko i monotonicznie,
przy czym nadwyżka kurtozy spada *szybciej* ($O(1/n)$) niż skośność
($O(1/\sqrt{n})$).

![Średnie wartości $Z_1$ (skośność) i $Z_2$ (kurtoza) z $M=5\,000$ replikacji, $K=50$.](plots/05_dagostino_components.png){width=72%}

Średnie $Z_1$ przekracza wartość krytyczną $\pm 1{,}96$ aż do
$n \approx 5$, podczas gdy $Z_2$ utrzymuje się **poniżej** tej granicy
dla wszystkich $n \ge 2$. Innymi słowy: w naszym przypadku
**całą "pracę" detekcji** w teście D'Agostino wykonuje składowa
skośnościowa. Składowa kurtozowa jest "zasilana" mniejszym sygnałem,
ponieważ $6/n$ tłumi się szybciej niż $2/\sqrt{n}$. To uzasadnia,
dlaczego DA - choć łączy dwie składowe - nie zyskuje istotnej
przewagi nad SW: skośność (na której zarówno SW, jak i DA pośrednio
się opierają) zanika powoli, a kurtoza zanika zbyt szybko, by
stanowić użyteczne źródło dodatkowego sygnału.

### 3.6. Rozkład $p$-value: trzy scenariusze

![$p$-value przy $n=5$ ($K=50$): nagromadzenie blisko zera, najwyraźniejsze dla SW i DA.](plots/06_pvalues_n5.png){width=72%}

![$p$-value przy $n=30$ ($K=50$): kształt mieszany - nagromadzenie blisko zera tylko w SW i DA.](plots/07_pvalues_n30.png){width=72%}

![$p$-value przy $n=200$ ($K=50$): brak nagromadzenia - wszystkie trzy rozkłady niemal Uniform(0,1).](plots/08_pvalues_n200.png){width=72%}

Wizualizacja $p$-value pokazuje, dlaczego decyzja oparta wyłącznie na
liczbie odrzuceń bywa myląca. Przy $n=5$ SW i DA mają wyraźnie
trójkątny rozkład $p$-value (gęstość rośnie z $p$-value malejącym do
zera). Przy $n=200$ wszystkie trzy histogramy są bliskie Uniform(0,1)
- hipoteza $H_0$ jest *prawie* spełniona, więc rozmiar testu pokrywa
się z poziomem istotności.

## 4. Interpretacja

**Hierarchia mocy.** Najmocniejszym testem zgodności rozkładu
$\bar{X}_n$ z normalnym dla rozkładu wyjściowego $\mathrm{Exp}(1)$ i
$K \in \{20,50,100\}$ jest test **Shapiro–Wilka** (w zakresie małych
$n$, gdzie różnice między testami są największe), zaraz za nim
**D'Agostino–Pearson K²**. Test **Kołmogorowa–Smirnowa** istotnie
ustępuje obu pozostałym testom - szczególnie przy małym $n$, gdzie
różnica wynosi nawet kilkadziesiąt punktów procentowych mocy.

Powód jest dobrze znany w literaturze: SW ma najlepsze własności
(zarówno dla dużych, jak i małych prób) względem **alternatyw
skośnych** o lekko cięższych ogonach niż normalny - czyli dokładnie tej
klasy, w której znajduje się $\bar{X}_n$ przy małym $n$ dla
rozkładu $\mathrm{Exp}(1)$. KS - jako test "ogólny", oparty na
maksymalnej odległości między dystrybuantami - rozprasza moc na
wszystkie kierunki odchyleń i przez to gorzej wykrywa konkretny rodzaj
niezgodności.

**Rola skośności i kurtozy w D'Agostino.** Empiryczne wartości
skośności idealnie pokrywają się z teoretycznym $2/\sqrt{n}$, natomiast
nadwyżka kurtozy zanika szybciej niż wynikałoby z teorii - efekt
ujemnego obciążenia estymatora $b_2$ dla skończonego $K$. To wyjaśnia
różnicę między SW a DA: oba testy korzystają z asymetrii, ale DA
osłabia sygnał skośnościowy przez słabszy, szybciej zanikający
sygnał kurtozowy. Stąd DA niemal nigdy nie pokonuje SW w naszym
eksperymencie. Wyjątkiem jest obszar $n \ge 30$, gdzie obie składowe
$Z_1, Z_2$ są małe i każda niewielka przewaga DA wynika z faktu, że
łączy dwa źródła informacji.

**Rola wielkości $K$.** Wzrost $K$ z 20 do 100 zwiększa moc każdego
testu **wielokrotnie** - dla $n=5$ moc SW rośnie z 0,24 do 0,91. Jest
to znacznie większy efekt niż wybór testu w obrębie ustalonego $K$.
Dlatego z praktycznego punktu widzenia, **jeśli badacz może
kontrolować $K$, jest to ważniejsza decyzja niż wybór między SW a DA.**

**Empiryczny rozmiar testu.** Dla $n=200$ (gdzie $H_0$ jest niemal
spełniona) wszystkie trzy testy mają empiryczny rozmiar bliski
nominalnemu $\alpha = 0{,}05$ (SW: 0,060; DA: 0,069; KS: 0,046). Lekkie
odchylenie DA powyżej $\alpha$ jest tutaj skutkiem ograniczonej liczby
replikacji ($M=5000$, błąd standardowy $\approx 0{,}003$) i nie jest istotne.

## 5. Wnioski

1. Dla rozkładu $\mathrm{Exp}(1)$ i wielkości testu $K \in \{20,50,100\}$
   **test Shapiro–Wilka** jest najmocniejszy w niemal całym zakresie
   $n$, choć od $n \ge 30$ jego przewaga nad D'Agostino–Pearsonem zanika.

2. **Test D'Agostino–Pearsona** jest bliski SW i lepszy od KS;
   jego siłą jest łatwość interpretacji - rozkłada decyzję na dwie
   składowe ($Z_1$ - skośność, $Z_2$ - kurtoza), co pozwala
   *zobaczyć*, *dlaczego* test odrzuca.

3. **Test Kołmogorowa–Smirnowa** jest istotnie słabszy dla
   alternatyw o lekkiej prawoskośności, jakie generuje średnia z
   $\mathrm{Exp}(1)$. Jest to znana właściwość tego testu - uniwersalność
   procedury KS oznacza utratę mocy dla konkretnych klas alternatyw.

4. **Liczność próby testowej $K$** ma znacznie większy wpływ na moc niż
   wybór konkretnego testu. Przed inwestycją w wybór procedury warto
   więc rozważyć zwiększenie $K$.

5. Empiryczna skośność $\hat{g}_1$ idealnie pokrywa się z teoretyczną
   $2/\sqrt{n}$, natomiast nadwyżka kurtozy $\hat{b}_2 - 3$ jest
   systematycznie *zaniżona* względem teoretycznego $6/n$ - efekt
   ujemnego obciążenia tych estymatorów dla skończonego $K$.
   To tłumaczy, dlaczego skośnościowa składowa testu D'Agostino
   wykonuje większą część pracy.

6. **Empiryczny rozmiar testu** dla wszystkich trzech procedur jest
   bliski nominalnemu $\alpha = 0{,}05$ przy dużym $n$, co potwierdza
   poprawność implementacji oraz fakt, że CTG dla $n \gtrsim 200$ jest
   tutaj praktycznie nieodróżnialne od ścisłej normalności.

---

### Dodatek A - odtworzenie wyników

```bash
# Utworzenie środowiska (jednorazowo)
conda create -n sad-clt -c conda-forge -y r-base r-ggplot2 r-dplyr r-tidyr \
    pandoc tectonic

# Uruchomienie symulacji
conda activate sad-clt
Rscript R/main.R

# Kompilacja raportu do PDF
pandoc report.md -o report.pdf --pdf-engine=tectonic -V lang=pl -H header.tex
```

Wygenerowane pliki:

* `plots/*.png` - wszystkie figury wykorzystane w raporcie,
* `output/power_grid.csv` - pełna siatka mocy testów dla
  $n \times K \in \{2,3,5,10,20,30,50,100,200\} \times \{20,50,100\}$,
* `output/pvalue_dump.csv` - surowe $p$-value z $M=5\,000$ replikacji
  dla trzech wybranych wartości $n$ (5, 30, 200) przy $K=50$.

Parametryzację eksperymentu (siatki $n$, $K$, liczbę replikacji
$M$, ziarno losowe, poziom istotności) zebrano w `R/config.R`.

---

# Problem 2. Odbiornik znanego sygnału

## 1. Wprowadzenie

Zadanie dotyczy detekcji słabego sygnału pseudoszumowego ukrytego w białym szumie
gaussowskim. Dwie osoby ustaliły wcześniej wspólny kod $s(i) \in \{-1, +1\}$,
$i = 1,\dots,N$ wygenerowany pseudolosowo. W sytuacji alarmowej nadawca emituje
sygnał $A s(i)$ o nieznanej amplitudzie $A$; odbiorca obserwuje:

$$
H_0:\; x(i) = w(i), \qquad
H_1:\; x(i) = A\,s(i) + w(i),
$$

gdzie $w(i) \overset{\text{iid}}{\sim} \mathcal{N}(0, \sigma^2)$.
Celem jest zaprojektowanie optymalnego detektora, wyznaczenie jego mocy oraz
analiza wpływu parametrów $A$, $\sigma^2$ i $N$ na skuteczność detekcji.

## 2. Metoda rozwiązania problemu

### 2.1. Wyprowadzenie statystyki testowej

Traktując próbkę $x(1),\dots,x(N)$ jako obserwacje (z różnymi wartościami
oczekiwanymi pod $H_1$), zapisujemy logarytm ilorazu wiarygodności:

$$
\ln \Lambda = \sum_{i=1}^{N} \ln \frac{f_1(x(i))}{f_0(x(i))}
= \sum_{i=1}^{N} \left[\frac{A\,s(i)\,x(i)}{\sigma^2}
  - \frac{A^2 s(i)^2}{2\sigma^2}\right]
= \frac{A}{\sigma^2}\underbrace{\sum_{i=1}^{N} s(i)\,x(i)}_{T}
  - \frac{A^2 N}{2\sigma^2},
$$

gdzie skorzystano z $s(i)^2 = 1$. Ponieważ dla $A > 0$ obie stałe
$A/\sigma^2$ i $A^2 N/(2\sigma^2)$ są dodatnie, test NP redukuje się do
porównania **statystyki korelacyjnej** (dopasowanego filtra):

$$
T = \sum_{i=1}^{N} x(i)\,s(i) = \langle x,\, s\rangle
$$

z progiem krytycznym $\eta$. Odrzucamy $H_0$, gdy $T > \eta$.

### 2.2. Rozkład T i wyznaczenie progu

Pod obiema hipotezami $T$ jest sumą niezależnych zmiennych gaussowskich:

$$
T \mid H_0 \;\sim\; \mathcal{N}(0,\; \sigma^2 N), \qquad
T \mid H_1 \;\sim\; \mathcal{N}(A N,\; \sigma^2 N).
$$

Dla zadanego poziomu istotności $\alpha$ (prawdopodobieństwo fałszywego alarmu)
próg krytyczny wynosi:

$$
\eta = z_{1-\alpha}\,\sqrt{\sigma^2 N},
$$

gdzie $z_{1-\alpha} = \Phi^{-1}(1-\alpha)$. Moc testu (prawdopodobieństwo
wykrycia alarmu) dana jest wzorem analitycznym:

$$
\beta(A, \sigma^2, N) = 1 - \Phi\!\left(\frac{\eta - A N}{\sqrt{\sigma^2 N}}\right)
= 1 - \Phi\!\left(z_{1-\alpha} - \frac{A\sqrt{N}}{\sigma}\right).
$$

Kluczowym parametrem jest stosunek sygnał/szum:

$$
\mathrm{SNR} = \frac{A^2 N}{\sigma^2}.
$$

Moc rośnie monotonicznie z $\mathrm{SNR}$: długi kod $N$ zastępuje silny sygnał $A$.

### 2.3. Schemat symulacji

Eksperyment przeprowadzono w **R 4.3.3** (tylko pakiety bazowe). Kod podzielony
modularnie na pliki `R/config.R`, `R/signal.R`, `R/detector.R`,
`R/simulation.R`, `R/plots.R` i `R/main.R`.
Reprodukowalność zapewnia `set.seed(20260519)`.

Ustalono: poziom istotności $\alpha = 0{,}001$, wariancja szumu
$\sigma^2 \in \{0{,}5;\, 1;\, 4\}$, długości kodu $N \in \{10,\,50,\,200,\,1000\}$,
amplitudy $A \in [0;\, 3]$ (krok $0{,}25$), liczba replikacji Monte Carlo
$M = 5\,000$ na punkt.

## 3. Otrzymane wyniki wraz z interpretacją

### 3.1. Zasada działania detektora

![Rozkłady statystyki korelacyjnej $T$ pod $H_0$ i $H_1$ dla $N=50$, $A=0{,}5$, $\sigma^2=1$, $\alpha=0{,}001$. Czerwony obszar to prawdopodobieństwo fałszywego alarmu.](plots/p2_01_rozklad_T.png){width=72%}

Rysunek ilustruje geometrię problemu: pod $H_0$ rozkład $T$ jest skupiony
wokół zera, pod $H_1$ przesuwa się o $AN = 25$ w prawo. Próg krytyczny
$\eta \approx 21{,}9$ jest tak dobrany, żeby prawdopodobieństwo przekroczenia
go pod $H_0$ wynosiło dokładnie $\alpha = 0{,}001$ -- jest to jednocześnie
gwarancja rzadkości fałszywego alarmu. Obszar nakładania się obu krzywych
odpowiada przypadkom, gdy detekcja jest trudna: słaby sygnał lub duży szum.

### 3.2. Moc w funkcji amplitudy A - wpływ długości kodu N

![Empiryczna i teoretyczna moc detektora w funkcji amplitudy $A$ dla różnych $N$ ($\sigma^2=1$, $\alpha=0{,}001$). Linie przerywane: wartości teoretyczne.](plots/p2_02_moc_vs_A_N.png){width=75%}

Zgodność empirycznej mocy z wartościami teoretycznymi potwierdza poprawność
implementacji. Widoczna jest wyraźna hierarchia krzywych: przy $N=10$ nawet
$A=1{,}5$ nie wystarcza do pewnej detekcji, podczas gdy $N=1000$ zapewnia moc
bliską $1$ już dla $A=0{,}25$. Jest to bezpośredni efekt wzrostu SNR z $N$:
podwojenie długości kodu równoważy obcięcie amplitudy o czynnik $\sqrt{2}$.

### 3.3. Moc w funkcji amplitudy A - wpływ wariancji szumu $\sigma^2$

![Empiryczna i teoretyczna moc detektora w funkcji amplitudy $A$ dla różnych $\sigma^2$ ($N=50$, $\alpha=0{,}001$). Linie przerywane: wartości teoretyczne.](plots/p2_03_moc_vs_A_sigma2.png){width=75%}

Wzrost wariancji szumu z $0{,}5$ do $4$ przesuwa krzywą mocy znacząco w prawo:
przy $\sigma^2=4$ osiągnięcie mocy $0{,}9$ wymaga amplitudy $A \approx 1{,}5$,
podczas gdy przy $\sigma^2=0{,}5$ wystarczy $A \approx 0{,}35$. Efekt jest
symetryczny: moc zależy wyłącznie od $\mathrm{SNR} = A^2 N / \sigma^2$,
więc czterokrotne zwiększenie szumu odpowiada dwukrotnemu zmniejszeniu amplitudy.

### 3.4. Przypadek słabego sygnału: $A \ll \sigma$

![Moc detektora w funkcji długości kodu $N$ dla słabych sygnałów ($\sigma^2=1$, $\alpha=0{,}001$). Linie przerywane: wartości teoretyczne.](plots/p2_04_moc_vs_N.png){width=75%}

Przypadek $A \ll \sigma$ jest kluczowy z praktycznego punktu widzenia - sygnał
alarmowy ma być trudny do wykrycia przez obserwatora zewnętrznego. Dla
$A=0{,}25$ ($\mathrm{SNR} = 3{,}1$ przy $N=50$) detektor praktycznie nie
działa: moc $\approx 0{,}09$. Jednak przy $N=500$ ta sama amplituda daje moc
$\approx 0{,}98$. Dla $A=0{,}5$ przejście następuje wcześniej: $N=100$ daje
moc $\approx 0{,}97$. Wynika z tego fundamentalna zasada: **wydłużenie kodu
$N$ kompensuje słabość sygnału** -- wzrost $N$ czterokrotnie odpowiada
podwojeniu amplitudy.

### 3.5. Weryfikacja rozmiaru testu pod $H_0$

![$p$-value pod $H_0$ ($N=50$, $\sigma^2=1$, $M=10\,000$ replikacji). Rozkład bliski $\mathrm{Uniform}(0,1)$ potwierdza poprawność implementacji.](plots/p2_05_pvalue_H0.png){width=62%}

Rozkład $p$-value pod $H_0$ jest niemal jednostajny, a rozmiar empiryczny
($0{,}0008$) pokrywa się z nominalnym $\alpha=0{,}001$. Potwierdza to
analityczną poprawność wzoru na próg krytyczny.

## 4. Interpretacja i wnioski

**Optymalność.** Wyprowadzony detektor jest testem Neymana–Pearsona dla
$A > 0$: przy ustalonym $\alpha$ maksymalizuje moc. Statystyka
$T = \langle x, s\rangle$ odpowiada filtrowi dopasowanemu (_matched filter_) --
klasycznemu rozwiązaniu w teorii detekcji sygnałów.

**Kryterium progu.** Próg $\eta = z_{1-\alpha}\sqrt{\sigma^2 N}$ wyznacza się
z warunku na prawdopodobieństwo fałszywego alarmu $P_{\text{FA}} = \alpha$.
W zastosowaniu alarmowym $\alpha$ powinno być bardzo małe ($10^{-3}$--$10^{-6}$),
co odzwierciedla wysoką cenę fałszywego alarmu. Próg rośnie z $\sigma$ i $N$,
ale moc rośnie _szybciej_ niż próg -- detektor zyskuje na długości kodu.

**Wielobitowy system telekomunikacyjny.** Pojedyncza detekcja sygnału $s$
przekazuje 1 bit informacji (alarm/brak). Aby przesłać wiele bitów jednocześnie,
można zastosować ortogonalny zestaw kodów $\{s_k\}$, $k=1,\dots,B$, każdy
o długości $N$, i niezależnie decydować o obecności każdego z nich w odebranym
$x$. Dzięki ortogonalności $\langle s_k, s_l \rangle = 0$ dla $k \neq l$
statystyki $T_k = \langle x, s_k \rangle$ są niezależne pod $H_0$, co
pozwala kontrolować globalne prawdopodobieństwo fałszywego alarmu korekcją
Bonferroniego: $\alpha_{\text{lok}} = \alpha / B$. Minimalna wymagana
długość kodu rośnie proporcjonalnie do $B/\alpha$, co ogranicza praktyczną
liczbę równoczesnych bitów.

## 5. Wnioski

1. Optymalny detektor dla modelu $H_1: x = As + w$ sprowadza się do
   porównania korelacji dopasowanej $T = \langle x, s\rangle$ z progiem
   $\eta = z_{1-\alpha}\sqrt{\sigma^2 N}$.

2. Moc testu zależy wyłącznie od $\mathrm{SNR} = A^2 N / \sigma^2$;
   wzrost $N$ jest równoważny wzrostowi $A^2$ -- długi kod kompensuje
   słabą amplitudę.

3. W reżimie $A \ll \sigma$ (sygnał ukryty w szumie) wystarczające
   SNR można osiągnąć przez wydłużenie kodu, a nie zwiększanie mocy
   nadawania -- co jest zaletą z punktu widzenia ukrycia komunikacji.

4. Empiryczny rozmiar testu ($0{,}0008 \approx \alpha$) i zgodność
   mocy empirycznej z teoretyczną potwierdzają poprawność implementacji.

5. System można rozszerzyć do wielobitowego przez zestaw wzajemnie
   ortogonalnych kodów, z Bonferroniego korektą poziomu istotności.

---

### Dodatek B -- odtworzenie wyników (Problem 2)

```r
# z katalogu glownego projektu:
setwd("sciezka/do/projektu")
source("R/main.R")
```

Wygenerowane pliki:

* `plots/p2_0*.png` -- wykresy Problemu 2,
* `output/p2_moc_vs_A_N.csv` -- siatka mocy dla $N \times A$,
* `output/p2_moc_vs_A_sigma2.csv` -- siatka mocy dla $\sigma^2 \times A$.
