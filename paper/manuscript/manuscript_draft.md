# Manuscript draft: release shape, not more grades, aligns single-application controlled-release urea with cereal N demand

Oct 2, 2026 · @Tianlin

## Title page

**Title.** Release shape, not more grades, aligns single-application controlled-release urea with nitrogen demand of rice, wheat and maize worldwide

**Authors.** **\[TO FILL: author names; corresponding author with e-mail; ORCID iDs\]**

**Affiliations.** **\[TO FILL: company R&D centre; academic partner institution(s)\]**

**Highlights**

- A model links static-water release parameters of coated urea to field N supply and crop demand.
- Urea plus one linear grade saves 12–28% N against equally safe split urea.
- Adding a second linear grade saves almost no further N (median 0%).
- A moderately S-shaped release (Weibull β ≈ 2–2.5) minimises N need and yield risk.
- Strongly lagged release (β > 3.5) raises shortfall risk in warm or cool years.

**Keywords.** controlled-release urea; polymer-coated urea; nitrogen synchrony; release kinetics; Weibull model; single basal application

*Note to authors (delete before submission).* Gaps are marked **\[TO FILL: …\]** (needs your data or a decision) and **\[VERIFY: …\]** (check against the source or rerun the model). In the Word version these are highlighted in yellow. Every number comes from model runs with literature parameters and must be recomputed after your products are calibrated (Section 2.8).

## Abstract

**Context.** A single basal application of controlled-release urea (CRU) blended with urea saves labour and, on average, raises yield and lowers nitrogen (N) losses. Yet products are still chosen by their labelled release period, and no framework links a product's release parameters to crop N demand across climates.

**Objective.** To derive, for rice, wheat and maize worldwide, the control period, urea share and release shape that meet crop N demand with the least N, and to test whether adding more product grades or changing the release shape matters more.

**Methods.** Release in 25 °C static water was described by a Weibull function whose shape parameter β separates near-linear (β ≈ 1.3) from sigmoidal (β ≥ 2) products. It was converted to field release with Arrhenius temperature scaling, integration over the diurnal cycle, a freezing cut-off and soil-moisture, placement and texture factors. Crop N demand followed growing-degree-day logistic curves calibrated at two growth stages. A daily mass balance of fertiliser-derived soil N was combined with a linear program that minimised total N subject to no shortfall in normal, 1.5 °C cooler and 1.5 °C warmer seasons. The framework was applied to 73 rice seasons in 50 zones, 33 wheat zones and 28 maize zones, with β scanned from 1.0 to 5.0. Environmental effects were estimated with IPCC emission factors and meta-analytic CRU response ratios.

**Results.** Urea plus one linear grade required a median 14% (rice), 12% (wheat) and 28% (maize) less N than an equally safe split-urea programme. A second linear grade saved a median 0% more (maximum 1.8%). Replacing the linear grade with one sigmoidal grade (β = 2.5) saved a further 11% (rice), 6% (wheat) and 7% (maize). At 85% of the N rate, the mean yield-loss proxy fell from 2.4% to 0.5% in rice, from 3.1% to 1.6% in wheat on freezing soils and from 3.8% to 2.2% in temperate maize as β rose from 1.3 to its optimum. The optimum was β ≈ 2.5 for rice, 2.0–2.5 for wheat and ≈ 2.0 for maize. Above β ≈ 3.5, worst-case shortfall exceeded that of current linear products.

**Conclusions.** Release shape, not the number of grades, is the main lever for matching single-application N supply to cereal demand. A moderately sigmoidal coating (β ≈ 2–2.5), with the control period tuned to climate, serves all three cereals.

**Implications.** The results give coating developers a quantitative design target. They also provide a validation protocol for field release, pairing soil-driven and weather-driven runs with leave-one-site-out calibration. **\[TO FILL: one sentence with field-validation results once burial-bag data exist.\]**

## 1. Introduction

Cereals receive most of the world's synthetic nitrogen, yet a large share of it is lost to air and water. These losses impose health, climate and ecosystem costs that often exceed the farm value of the fertiliser (Zhang et al., 2015; Gu et al., 2023). Much of the loss arises because N is applied when the crop takes up little of it. Basal urea dissolves within days, while uptake peaks weeks later, between stem elongation and flowering. Splitting applications improves synchrony but costs labour, and labour is increasingly scarce in the major cereal regions. County-scale analyses in China suggest that better-timed and localised N management could halve fertiliser input without yield loss (Liu et al., 2024).

Controlled-release urea offers synchrony without extra field operations. A global meta-analysis of 240 studies found that replacing readily available N with controlled-release N raised cereal yields by about 5% and cut N losses by 33–49% (Zhang et al., 2024). Similar gains were reported for rice (Jiang et al., 2021). Single basal applications of CRU, often blended with urea to cover early demand and lower cost, matched or exceeded split urea in rice and wheat (Chen et al., 2020; Zhang et al., 2022). The same meta-analyses show that the benefit varies widely with N rate, the CRU share of the blend, climate and coating type (Zhang et al., 2024). This variation suggests that the choice of product, not the technology itself, often limits performance.

Release from polymer-coated urea is diffusion-controlled. It follows a lag, a near-linear phase and a decay phase whose relative lengths depend on coating permeability and granule geometry (Shaviv et al., 2003; Lipin et al., 2023). Temperature controls the release rate through an Arrhenius-type dependence (Gandeza et al., 1991; Du et al., 2006; Adams et al., 2013). Field release departs from the static-water label under fluctuating temperature and surface placement (Ransom et al., 2020). Manufacturers nevertheless describe products by a single number, the days to reach 80% release in water at 25 °C (GB/T 23348; **\[VERIFY: equivalent ISO/EN standard, e.g. EN 13266\]**). Agronomists choose among these labels largely by trial and error.

Two questions about product design remain open. First, should a blend contain several release periods, as is common in practice, or is one enough? Second, what release shape best matches crop demand? Sigmoidal products exist (Aoki, 2019), and in South China a sigmoidal-CRU blend tracked rice uptake better than a parabolic one (Huang et al., 2019). However, no study has asked which shape is optimal, or whether a better shape outperforms more grades. Demand-informed CRU design has been modelled for sugarcane at five Australian sites (Zhao et al., 2017), and crop models have optimised rate and release period for single sites (Ren et al., 2024). No framework yet connects release parameters a coating engineer can tune to crop demand across the climates where cereals are grown.

Here we develop such a framework and apply it to rice, wheat and maize worldwide. Our objectives were to: (i) build a transparent model from static-water release parameters to field N supply and crop demand; (ii) derive the optimal control period, urea share and number of grades for each production system; (iii) identify the release shape (Weibull β) that minimises N need and yield risk, including in abnormally warm or cool seasons; and (iv) propose a field protocol to validate and calibrate the release model. We hypothesised that release shape, not the number of grades, is the main determinant of how closely a single application can follow crop demand.

## 2. Materials and methods

### 2.1. Release in static water

Cumulative release of a coated-urea grade in water at 25 °C was described by a Weibull function with an initial burst:

```latex
F(\tau) = f_0' + (1-f_0')\left[1-\exp\!\left(-(\tau/\lambda)^{\beta}\right)\right]
```

where τ is time (d) at 25 °C and f₀′ is the fraction released within the first day. The control period D follows the Chinese standard definition: the number of days to 80% release in 25 °C water (GB/T 23348). The scale parameter λ therefore follows from D, f₀′ and β:

```latex
\lambda = D\,\Big/\left[\ln\frac{1-f_0'}{0.2}\right]^{1/\beta}
```

The shape parameter β distinguishes near-linear release (β ≈ 1–1.5; a short lag, then near-constant rate) from sigmoidal release (β ≥ 2; a distinct lag followed by rapid release). The default for current products was β = 1.3 and f₀ = 3% **\[TO FILL: replace with values fitted to your grades, Section 2.8\]**. Granules damaged in handling and blending (2%) were assumed to dissolve immediately: f₀′ = f₀ + (1 − f₀) × 0.02. Product grades were D = 30–360 d in 10-d steps (34 grades), matching the range manufactured by the authors' company **\[VERIFY: confirm product range\]**.

### 2.2. Field release

Field release was obtained by accumulating 25 °C-equivalent days at the coating:

```latex
\frac{d\tau}{dt} = f_{soil}\; a_P\; a_S\; \left[1 - s\,x(\psi)\right]\; \overline{\exp\!\left[\frac{E_a}{R}\left(\frac{1}{298.15}-\frac{1}{T+273.15}\right)\right] r(T)}
```

The overbar denotes the mean over a sinusoidal diurnal cycle of the coating temperature T. Because the Arrhenius term is convex, a wider daily range accelerates release at the same mean temperature, consistent with faster release under fluctuating than constant temperature (Ransom et al., 2020). The freezing term r(T) falls linearly from 1 at 3 °C to 0 at 0 °C, so release stops in frozen soil. The apparent activation energy Ea was 50 kJ mol⁻¹ (Q₁₀ ≈ 2 at 20–30 °C; Gandeza et al., 1991; Du et al., 2006 reported 37–46 kJ mol⁻¹).

The moisture factor reduces release as soil dries: x(ψ) rises linearly in log₁₀|ψ| from 0 at −33 kPa to 1 at −1500 kPa, and the sensitivity s = 0.25 delays a 60-d product by about 27 d at −1000 kPa and 20 °C **\[VERIFY: against Verburg et al. 2020 and add reference\]**. Placement factors a\_P were 1.0 (incorporated), 0.85 (band or deep placement) and 2.5 (unincorporated surface application, with +2 °C and a ±10 °C diurnal range; Ransom et al., 2020). Texture factors a\_S were 1.15 (clay), 1.0 (loam) and 0.95 (sand) **\[VERIFY: add Golden et al. 2011 reference\]**. The soil factor f\_soil = 0.90 accounts for slower release in moist soil than in free water.

Coating temperature and moisture regime came from cropping-system presets. Coating temperature was air temperature + 1 °C; under air temperatures below 0 °C it was 0.4 × air temperature + 0.5 °C to represent snow and soil buffering. Film mulch added 3.5 °C for 60 d, decreasing to 1 °C over the next 40 d. Moisture regimes were flooded (ψ = 0), alternate wetting and drying, irrigated upland (−10 to −80 kPa), humid rainfed (−10 to −300 kPa) and semi-arid rainfed (−50 to −1200 kPa).

### 2.3. Crop N demand

Cumulative N uptake was a normalised logistic function of growing-degree-day (GDD) progress x = GDD(t)/GDD(maturity), fitted through two calibration points per crop (Table 1). The points were panicle initiation and heading for rice, jointing and anthesis for wheat, and V10 and silking for maize. For maize, about 60% of seasonal N is taken up by silking (Ciampitti and Vyn, 2012; Bender et al., 2013). Total uptake was target yield × N requirement per tonne of grain (rice 17–19, wheat 25–29, maize 19 kg N t⁻¹). The N to be supplied by fertiliser was total uptake minus indigenous soil supply (the N uptake of zero-N plots), with a floor of 12–15 kg N ha⁻¹ **\[VERIFY: references for N requirement per tonne and stage fractions of rice and wheat\]**.

### 2.4. Fertiliser-derived soil N pool

Fertiliser N released into the soil entered a pool P that lost a fraction k per day and supplied crop uptake with root capture efficiency η:

```latex
P(t+1) = P(t)\,(1-k) + R(t) - \frac{U(t)}{\eta}
```

R(t) is daily release from the blend and U(t) daily fertiliser-N uptake. Urea hydrolysed with a time constant of 2–3 d after an initial loss of 10–15%, mainly ammonia volatilisation. For rice k = 0.013 d⁻¹ (×1.15 for dry-seeded rice). For upland crops k = k₂₀ · 2^((T−20)/10) · w, with k₂₀ = 0.010 d⁻¹ and w > 1 in wet regions. η was 0.80. A season was without shortfall when P(t) never fell below the larger of: next 7 days of demand ÷ η, and an early-season reserve of 10–20% of fertiliser demand until panicle initiation, jointing or V10 **\[TO FILL: calibrate k and η against ¹⁵N or apparent-recovery data; current values give urea recoveries of about 40–50%\]**.

### 2.5. Blend optimisation

For each season, the contribution of 1 kg N of each product k to the pool, G\_k(t), was computed from the pool equation without uptake. The minimum-N blend then solved:

```latex
\min \sum_k c_k x_k \quad \text{s.t.} \quad \sum_k x_k\, G_k^{(\Delta T)}(t) \ge N^{(\Delta T)}(t)\quad \forall t,\ \Delta T \in \{-1.5,\,0,\,+1.5\}
```

Here x\_k is the N rate of product k and N(t) is cumulative demand plus the no-shortfall floor. Cost weights were c = 1 for urea and 1.02 for CRU, so that urea is preferred when equally effective. Stacking constraints for a normal season and seasons 1.5 °C cooler and warmer makes one recipe robust to typical interannual variation. Candidate blends were urea plus one or two grades; grades releasing < 30% by maturity were excluded. Shares were rounded to 5%, by enumerating nearby integer combinations and keeping the one with the least total N. The total was rounded up to 5 kg N ha⁻¹. A single-grade blend was preferred if it needed ≤ 3% more N than the best two-grade blend. The benchmark was urea split into two or three doses (rice 40/30/30% at transplanting, tillering and panicle initiation; wheat 50/50% at sowing and jointing; maize 40/60% at sowing and V10), with the same no-shortfall criterion.

The rice, wheat and maize recipe maps (Sections 3.2–3.3) use crop-specific release modules. The rice map uses coating temperature = air + 1 °C without diurnal integration and is optimised for the normal season only. The wheat map adds the freezing cut-off with f\_soil = 0.80. The maize map uses the full field model of Section 2.2. **\[TO FILL: rerun rice and wheat maps with the unified field model before submission so that all results share one release model; differences are expected to be small.\]**

### 2.6. Release-shape scan

To isolate the effect of release shape, every season was re-optimised with the unified field model (Section 2.2) for β = 1.0, 1.3, 1.6, 2.0, 2.5, 3.0, 3.5, 4.0 and 5.0. Each β was free to choose its best D and urea share. Two comparisons were made. In the equal-yield comparison, the minimum N for no shortfall (Section 2.5) was expressed relative to β = 1.3. In the equal-N comparison, all β received the same N: 100% or 85% of the β = 1.3 requirement. Each β chose the D and urea share minimising mean unmet uptake over the three design seasons, and was then tested under five conditions: normal, 3 °C cooler, 3 °C warmer, field release 25% slower and 33% faster. Unmet uptake was obtained from the pool equation with uptake limited by pool size. The yield-loss proxy was unmet uptake ÷ total crop N uptake. This assumes yield is proportional to N uptake in the N-limited range and is therefore an upper-bound indicator rather than a yield prediction.

### 2.7. Environmental effects

NH₃ volatilisation, direct and indirect N₂O, and N leaching and runoff were estimated as N rate × emission factor × (1 − c × RR). Emission factors followed IPCC (2019) and were adjusted for climate and soil pH. c is the CRU share of the blend, and RR the meta-analytic reduction by CRU relative to urea at equal N (rice: Jiang et al., 2021; upland crops: Zhang et al., 2024; Lu et al., 2023 **\[VERIFY: full reference for Lu 2023\]**). Greenhouse gases included fertiliser manufacture (8 kg CO₂e kg⁻¹ N in China, 4 elsewhere), coating polymer (+0.43 kg CO₂e kg⁻¹ CRU-N), N₂O, and paddy CH₄ (IPCC 2019 daily factors). Damage costs used unit costs for NH₃ and leached N following Gu et al. (2023) **\[VERIFY: unit damage costs and source\]**, and a carbon price of 50 USD t⁻¹ CO₂e. Uncertainty was propagated by Monte Carlo sampling (2000 draws, triangular distributions), using common random numbers across scenarios. Because RR does not depend on release shape, differences among β arise only from N rate and CRU share.

### 2.8. Case data, calibration and validation

The case database comprised 50 rice zones (73 seasons), 33 wheat zones and 28 maize zones covering the main production regions (Fig. 2). For each zone it held a representative station with 12 monthly mean temperatures, the crop calendar, target yield, indigenous N supply and local N rate, drawn from FAO/GEOGLAM crop calendars, national statistics and the literature (Supplementary Table S3). Daily temperature was interpolated from monthly normals **\[TO FILL: replace with 20-year AgERA5 or NASA POWER daily data; report interannual variability instead of fixed ±1.5 °C\]**.

Product calibration: static-water release of **\[TO FILL: number\]** grades was measured at 15, 25 and 35 °C and fitted for D, β, f₀ and Ea by least squares (code in Supplementary S7) **\[TO FILL: results\]**.

Field validation used buried mesh bags at **\[TO FILL: number\]** sites (**\[TO FILL: sites, crops, seasons\]**), retrieved at 0–240 d, with N remaining measured by **\[TO FILL: method\]**. The model was run twice. Run A used measured hourly soil temperature at bag depth and measured soil moisture converted to matric potential. Run B used daily air temperature and rainfall, with system presets and a single-layer water balance. The difference between runs quantifies the cost of using only weather data. f\_soil, s, and the placement and texture factors were calibrated on run A with weak priors, using leave-one-site-out cross-validation. Acceptance criteria were RMSE ≤ 8 percentage points, error in the day of 80% release ≤ 10 d or 15%, and run-B RMSE ≤ 1.5 × run-A RMSE. The workflow was tested on simulated data with known parameters (Fig. S1) **\[TO FILL: replace with field results\]**.

All analyses used Python 3.11 with NumPy and SciPy (HiGHS linear-programming solver). Code and data are available at **\[TO FILL: Zenodo DOI\]**.

## 3. Results

### 3.1. From label to field: a worked example

In water at 25 °C, linear (β = 1.3) and sigmoidal (β = 2.5) grades with the same D = 90 d both reach 80% at day 90. Before that, the linear grade releases more: 54% versus 33% by day 50 (Fig. 1a). For single-season hybrid indica rice at Wuhan (transplanted in **\[VERIFY: date\]**, 120 d to maturity, 200 kg N ha⁻¹ local practice), the minimum-N linear blend was urea 20% + CR90 80% at 195 kg N ha⁻¹. The equally safe three-way urea split needed 215 kg N ha⁻¹. Field release of the linear grade ran far ahead of demand during the first 60 d (Fig. 1b), and the soil pool peaked at about 65 kg N ha⁻¹ around day 45, well above the no-shortfall floor (Fig. 1c). This early surplus is exposed to loss, so the total must be raised until the pool still covers demand at heading. A sigmoidal blend (urea 25% + S100 75%) held the pool close to the floor and met the same constraints with 165 kg N ha⁻¹, 15% less.

### 3.2. Global recipes for current linear products

Urea plus one linear grade satisfied the no-shortfall constraints in every season: all 73 rice seasons, 30 of 33 wheat zones and 26 of 28 maize zones used a single grade (Fig. 2; Table 2). Rice used CR40–CR130, mostly CR70–CR90, with 15–25% urea. Wheat used short grades (CR30–CR50) with 15–35% urea, because cool autumn and frozen winter soil slow release and over-winter losses are small. Maize used CR50–CR100 with 0–10% urea: about 5% urea was enough to cover seedling demand, and dropping it raised the N requirement by up to 37% (NE China: 175 vs. 240 kg N ha⁻¹).

Against equally safe split urea, the blends saved a median of 14.3% N in rice, 11.8% in wheat and 27.5% in maize (Fig. 3a). Compared with reported local rates, the recommended totals were a median 11% lower in rice, 4% higher in wheat and 2% lower in maize. Median modelled recovery of fertiliser N was 0.47 (rice), 0.53 (wheat) and 0.55 (maize). Median N still unreleased at maturity was 3–7% **\[VERIFY: recompute after harmonising release modules, Section 2.5\]**.

### 3.3. More grades or a different shape?

Allowing a second linear grade hardly changed the requirement: the median saving was 0.0% in both rice and wheat (maximum 0.2% and 1.8%; Fig. 3b). Replacing the linear grade with one sigmoidal grade (β = 2.5) saved a median 14.6% in rice and 5.7% in wheat in the shape comparison. Combining one linear and one sigmoidal grade added little to the sigmoidal grade alone (14.9% and 7.2%). In the maps, the sigmoidal single-grade scenario reduced N by a median 11.4% (rice), 5.9% (wheat) and 7.4% (maize; β = 2.0) **\[VERIFY: harmonise the two rice estimates (14.6% vs 11.4%), which differ in robustness constraints\]**.

### 3.4. Optimal release shape

The β scan gave a consistent optimum (Fig. 4). At equal yield, the N needed for no shortfall fell from β = 1.0 to a minimum at β = 2.5–3.0 in rice (87.8% of the β = 1.3 requirement) and wheat on freezing soils (90.1–90.9%). The minimum was at β = 2.0 in temperate (92.8%) and tropical maize (90.2%). Beyond the minimum the requirement rose again, exceeding that of β = 1.3 for temperate maize at β ≥ 3.5.

At equal N, the advantage of a sigmoidal release appeared when N was limiting. At 85% of the β = 1.3 requirement, the mean yield-loss proxy across the five test conditions fell from 2.4% to 0.5% in rice, from 3.1% to 1.6% in wheat on freezing soils, from 3.5% to 2.1% in other wheat, from 3.8% to 2.2% in temperate maize and from 2.9% to 1.3% in tropical maize (Fig. 4a). The optimum was β = 2.0–2.5 in every group. Worst-case losses showed the same pattern but rose steeply above β ≈ 3 (Fig. 4b). At β = 5, worst-case losses exceeded those of β = 1.3 in all five groups. At 100% N, shortfalls were rare for every β: the mean worst-case proxy was 0.9–1.6%, and slightly higher at the optimum β than at β = 1.3 in wheat (1.6% vs 1.0%) and maize (1.0% vs 0.9%).

Fig. 5 shows why strong lags carry risk. In NE China spring maize, the β = 1.3 recipe (urea 5% + D70, 172 kg N ha⁻¹) kept a large pool in all three seasons. The β = 2.0 recipe (urea 8% + D60, 159 kg N ha⁻¹) held the pool lower and ran 2.1% short in a 3 °C warmer season, against 0.4% for β = 1.3. A strong lag (β = 4.0) needed 37% urea to bridge the lag phase, more N in total (179 kg N ha⁻¹), and still ran 2.5% short when warm weather advanced release before silking.

### 3.5. Environmental effects

At equal yield, moving from β = 1.3 to the optimum reduced total N by a median 12.2% (rice), 6.8% (wheat) and 7.6% (maize; Fig. 6). NH₃ fell by 9.7%, 5.5% and 6.0%; N₂O by 10.5%, 5.4% and 6.4%; leaching and runoff by 11.2%, 5.1% and 6.1%; and damage costs by 10.4%, 5.6% and 6.3%. Greenhouse-gas emissions fell by 6–7% in wheat and maize but only 1.2% in rice, where paddy CH₄, independent of N rate in our approach, dominates. Because the response ratios do not depend on release shape, these changes reflect N saved, not a direct shape effect **\[TO FILL: if available, field measurements of NH₃/N₂O under sigmoidal vs linear CRU\]**.

## 4. Discussion

### 4.1. Why release shape matters more than the number of grades

Blends of several release periods are often promoted as a way to follow crop demand more closely. Our results suggest this rarely works with near-linear products. A second linear grade reduced the N requirement by a median 0% because grades of the same shape produce nested supply curves. Any two of them span almost the same set of attainable supply curves as one grade plus urea. In a linear program the optimum lies at a vertex of the feasible region, so an extra, nearly collinear option is seldom used. What limits a linear grade is its front-loaded release. Much N is released when demand is low and pool losses accumulate, and the total must be raised until enough remains at peak uptake (Fig. 1). A lag phase shifts release towards the peak-demand window, which no combination of linear grades can reproduce. This agrees with field observations in South China, where sigmoidal CRU in a urea blend tracked rice N uptake more closely than parabolic CRU and gave slightly higher yields in both rice seasons (Huang et al., 2019).

### 4.2. Why a moderate S-shape is optimal

The optimum at β ≈ 2–2.5 reflects a trade-off between synchrony and robustness. Stronger lags concentrate release in a shorter window. When temperature or field conditions shift that window, the error falls on the peak-demand period and cannot be compensated by the small urea share (Fig. 5). The 37% urea that a β = 4 recipe needed to bridge its lag phase also re-introduced the early losses that a sigmoidal release is meant to avoid. Maize favoured a lower β than rice. Two features of maize systems may explain this: spring warming already makes linear products release slowly at first and faster later, and maize takes up about 40% of its N after silking, so demand is spread over a longer period (Ciampitti and Vyn, 2012). Wheat on freezing soils benefited from a lag that keeps N in the granule over winter, but its optimum was capped by the risk of an early spring. **\[VERIFY: test these explanations with a sensitivity run that removes spring warming or post-silking uptake.\]**

The sigmoidal advantage was largest when N was limiting. At full N, every β rarely fell short, and the optimum β was slightly less robust than β = 1.3 under extreme conditions. Sigmoidal products should therefore be positioned as a means to cut N rates without yield penalty, rather than to raise yields at current rates. This agrees with meta-analytic evidence that CRU yield gains depend on N rate (Zhang et al., 2024).

### 4.3. Implications for coating design and product portfolios

The results point to a simple design target: a coating giving β ≈ 2–2.5 in 25 °C water, with D tuned by climate. In practical terms, a β = 2.5 product releases about 10% of its N by 0.3 D and 32% by 0.55 D (β = 2.0: 16% and 40%), compared with 30% and 53% for β = 1.3. These two checkpoints can serve as a quality-control specification for S-type grades. Diffusion models link the lag phase to the water permeability of the coating and the linear phase to solute permeability (Shaviv et al., 2003; Lipin et al., 2023). This gives coating chemists a physical handle on β independent of coating thickness, which mainly sets D. One sigmoidal coating platform, offered in a few control periods and blended with 5–25% urea, could serve rice, wheat and maize. A larger portfolio of linear grades would add complexity without saving N. **\[TO FILL: comment on cost of sigmoidal coatings relative to linear ones, from your production data.\]**

### 4.4. Limitations and next steps

The framework has five main limitations.

1. **Parameters.** All release, loss and uptake parameters are literature defaults. Product calibration and field validation (Section 2.8) are needed before the absolute N rates can be used for recommendations **\[TO FILL: report calibration and validation results here\]**.
2. **Climate input.** Monthly normals and fixed ±1.5 °C and ±3 °C offsets stand in for interannual variability; daily reanalysis data would allow season-by-season tests.
3. **Spatial representation.** Each zone is represented by one station and typical values of yield and indigenous N supply, so within-zone variation in soil, variety and management is not captured.
4. **Yield.** The yield-loss proxy is not a yield prediction; coupling with a crop model such as APSIM or DSSAT, as done for single sites by Ren et al. (2024), would allow yield to be simulated directly.
5. **Environmental response.** CRU response ratios do not depend on release shape, so environmental gains are driven by N saved. Measurements comparing sigmoidal and linear CRU at equal N are needed.

The relative conclusions are expected to be more robust than the absolute values, because they arise from the shapes of the supply and demand curves. These conclusions are that one grade suffices and that β ≈ 2–2.5 is optimal, with risk rising at higher β. **\[TO FILL: sensitivity analysis (Ea 38/65 kJ mol⁻¹, f\_soil ±20%, k ±50%, η 0.7/0.9) confirming that the optimal β stays within 2.0–2.5; Supplementary S5.\]**

## 5. Conclusions

A transparent model linking static-water release parameters to field N supply and crop demand was applied to rice, wheat and maize systems worldwide. Urea plus a single linear grade met demand with 12–28% less N than an equally safe split-urea programme, and adding a second linear grade gave almost no further saving. Release shape was the main lever. A moderately sigmoidal release (Weibull β ≈ 2.5 for rice, 2.0–2.5 for wheat and ≈ 2.0 for maize) reduced the N needed for no shortfall by 7–12%. At 15% less N, it roughly halved the yield-loss proxy. Strong lags (β > 3.5) increased risk in abnormal seasons. Coating development should therefore target a moderate S-shape rather than a wider range of linear grades, and field validation of release models should pair soil-driven and weather-driven runs. **\[TO FILL: one sentence on field-validation outcome.\]**

## Declarations

**CRediT authorship contribution statement.** **\[TO FILL: e.g. A: Conceptualization, Resources, Funding acquisition; B: Methodology, Software, Formal analysis, Visualization, Writing – original draft; C: Investigation (laboratory and field validation), Data curation; D: Supervision, Writing – review & editing.\]**

**Declaration of competing interest.** **\[TO FILL: e.g. “A, B and C are employed by \[company\], which manufactures controlled-release urea. The company had no role in the interpretation of the model results. D declares no competing interests.”\]**

**Funding.** **\[TO FILL: grant numbers or “This work was funded by \[company\].”\]**

**Declaration of generative AI and AI-assisted technologies in the writing process.** During the preparation of this work the authors used Claude (Anthropic) to develop model code, generate figures and draft the manuscript. After using this tool, the authors reviewed and edited the content as needed and take full responsibility for the content of the publication. **\[VERIFY: wording against the current Elsevier policy.\]**

**Data and code availability.** Model code, the zone database, all per-season results and the scripts that produce every figure are archived at **\[TO FILL: Zenodo DOI\]** and maintained at **\[TO FILL: public repository URL\]**. Interactive recipe maps for rice, wheat and maize are provided as supplementary web resources.

**Acknowledgements.** **\[TO FILL\]**

## Tables

**Table 1.** Model parameters, default values and sources. Values for rice / wheat / maize where they differ.

| Parameter | Symbol | Default | Source |
| --- | --- | --- | --- |
| Control period (80% release, 25 °C water) | D | 30–360 d, 10-d steps | GB/T 23348 |
| Weibull shape, current products | β | 1.3 **\[TO FILL: fitted\]** | this study |
| Initial burst | f₀ | 3% **\[TO FILL: fitted\]** | this study |
| Granule damage | — | 2% | assumed |
| Apparent activation energy | Ea | 50 kJ mol⁻¹ (sensitivity 38, 65) | Gandeza et al., 1991; Du et al., 2006 |
| Soil release factor | f\_soil | 0.90 / 0.80 / 0.90 | assumed **\[TO FILL: calibrated\]** |
| Moisture sensitivity | s | 0.25 | **\[VERIFY: Verburg et al., 2020\]** |
| Placement factor (incorporated / band / surface) | a\_P | 1.0 / 0.85 / 2.5 | Ransom et al., 2020 |
| Texture factor (clay / loam / sand) | a\_S | 1.15 / 1.0 / 0.95 | **\[VERIFY: Golden et al., 2011\]** |
| Coating temperature offset | — | +1 °C; film +3.5 °C for 60 d | assumed |
| Urea hydrolysis time constant | — | 2 / 3 / 3 d | assumed |
| Initial urea loss | — | 15% / 10% / 10% | **\[VERIFY: source\]** |
| Daily pool loss rate (20 °C for upland) | k | 0.013 / 0.010 / 0.010 d⁻¹ | assumed **\[TO FILL: calibrated\]** |
| Root capture efficiency | η | 0.80 | assumed |
| No-shortfall buffer | — | 7 d of demand | this study |
| Early-season reserve | — | 20% / 15% / 10% of fertiliser demand | this study |
| N requirement per t grain | — | 17–19 / 25–29 / 19 kg N t⁻¹ | Ciampitti and Vyn, 2012 **\[VERIFY: rice, wheat sources\]** |
| Uptake at calibration stages | — | rice PI 35–40%, heading 78–85%; wheat jointing 25–32%, anthesis 82–88%; maize V10 25%, silking 60% | **\[VERIFY: sources\]** |
| GDD base / cap | — | 10 / 30 °C (rice, maize); 0 / 30 °C (wheat) | — |
| Robust design band | ΔT | −1.5, 0, +1.5 °C | this study |

**Table 2.** Recommended blends by crop for current linear products, and the optimal release shape.

| Crop | Zones | Seasons | Grades used (d) | Urea share (%) | Median total N (kg ha⁻¹) | N saved vs. split urea (median %) | N saved by one S-type grade (median %) | β with lowest yield-loss proxy | β with lowest N need |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Rice | 50 | 73 | 40–130 | 15–25 | 105 | 14.3 | 11.4 | 2.5 | 2.5 |
| Wheat | 33 | 33 | 30–50 | 15–35 | 125 | 11.8 | 5.9 | 2.0 | 3.0 |
| Maize | 28 | 28 | 50–100 | 0–10 | 165 | 27.5 | 7.4 | 2.0 | 2.0 |

## Figure captions

Files: `paper/figures/` in the project repository (vector PDF and 600 dpi PNG).

**Fig. 1.** From label to field for single-season hybrid indica rice at Wuhan (normal season). (a) Cumulative release in 25 °C water of a linear (β = 1.3) and a sigmoidal (β = 2.5) grade with the same control period D = 90 d. (b) Cumulative field release of the minimum-N linear blend (urea 20% + CR90, 195 kg N ha⁻¹) and the sigmoidal blend (urea 25% + S100, 165 kg N ha⁻¹) against cumulative crop demand for fertiliser N. PI, panicle initiation. (c) Fertiliser-derived soil N pool for both blends against the no-shortfall floor.

**Fig. 2.** Recommended single-application blends for current linear products. Colour: control period of the main CRU grade. Symbol area: total N rate. (a) 73 rice seasons in 50 zones; (b) 33 wheat zones; (c) 28 maize zones. Where a zone has several seasons, symbols are offset slightly.

**Fig. 3.** (a) N saved by urea plus one linear grade relative to an equally safe split-urea programme. (b) Further N saved, relative to urea plus one linear grade, by a second linear grade, one sigmoidal grade (β = 2.5), or one linear plus one sigmoidal grade, for rice and wheat. Boxes: interquartile range and median; whiskers: 1.5 × IQR; points: seasons.

**Fig. 4.** Effect of release shape (Weibull β) for five crop groups. (a) Mean yield-loss proxy over five test conditions at 85% of the β = 1.3 N requirement. (b) Worst-case yield-loss proxy among those conditions. (c) Median N needed for no shortfall in normal and ±1.5 °C seasons, relative to β = 1.3. Dotted line: current products (β = 1.3); grey band: β = 2.0–2.5. The yield-loss proxy is unmet N uptake divided by total crop N uptake.

**Fig. 5.** Fertiliser-derived soil N pool in NE China spring maize for recipes designed with β = 1.3 (a), 2.0 (b) and 4.0 (c), each optimised to avoid shortfall within ±1.5 °C, in a normal season and seasons 3 °C cooler and warmer. Legends give the yield-loss proxy for each season.

**Fig. 6.** Change in total N, NH₃ volatilisation, N₂O emission, leaching and runoff, greenhouse-gas emission and damage cost when moving from β = 1.3 to the optimal β at equal yield (no shortfall). Bars: median across seasons; error bars: interquartile range.

**Fig. S1.** Validation workflow demonstrated on simulated burial-bag data from four synthetic sites with known parameters. These are not field measurements. Predicted versus observed cumulative release for run A (measured soil temperature and moisture) and run B (air temperature and rainfall), (a) with literature parameters and (b) after leave-one-site-out calibration. Dashed lines: ±8 percentage points. **\[TO FILL: replace with field data when available and move to the main text.\]**

## References

**\[VERIFY: complete author lists, volumes, pages and DOIs for every entry before submission; links point to the record used.\]**

- Adams, C., Frantz, J., Bugbee, B., 2013. [Macro- and micronutrient-release characteristics of three polymer-coated fertilizers: theory and measurements](https://consensus.app/papers/details/dbe6206341e154a69f9de2406a80c429/). J. Plant Nutr. Soil Sci.
- Aoki, H., 2019. [Development and application of coated fertilizer in Japan](https://consensus.app/papers/details/521a216a43f757de8ea7c65b6a552d49/). J. Food Sci. Eng.
- Bender, R.R., Haegele, J.W., Ruffo, M.L., Below, F.E., 2013. Nutrient uptake, partitioning, and remobilization in modern, transgenic insect-protected maize hybrids. Agron. J. 105, 161–170. **\[VERIFY\]**
- Chen, Z., et al., 2020. [Impact of controlled-release urea on rice yield, nitrogen use efficiency and soil fertility in a single rice cropping system](https://consensus.app/papers/details/604731fbb70e56b88dec27a77026fb05/). Sci. Rep.
- Ciampitti, I.A., Vyn, T.J., 2012. Physiological perspectives of changes over time in maize yield dependency on nitrogen uptake and associated nitrogen efficiencies: a review. Field Crops Res. 133, 48–67. **\[VERIFY\]**
- Du, C.-W., Zhou, J.-M., Shaviv, A., 2006. [Release characteristics of nutrients from polymer-coated compound controlled release fertilizers](https://consensus.app/papers/details/9196b693f08c58049344b53f26ecfe3c/). J. Polym. Environ.
- Gandeza, A.T., Shoji, S., Yamada, I., 1991. [Simulation of crop response to polyolefin-coated urea: I. Field dissolution](https://consensus.app/papers/details/4272983ef1f152748b2f6ab05577724b/). Soil Sci. Soc. Am. J.
- GB/T 23348, 2009. Slow release fertilizer. Standardization Administration of China. **\[VERIFY: current edition\]**
- Golden, B., et al., 2011. **\[TO FILL: full reference for texture effect on PCU release\]**
- Gu, B., et al., 2023. [Cost-effective mitigation of nitrogen pollution from global croplands](https://consensus.app/papers/details/8d2668d1febf5d088072708e99cb2101/). Nature.
- Huang, Q.-Y., et al., 2019. [Seasonal differences in N release dynamic of controlled-released urea in paddy field and its impact on the growth of rice under double rice cropping system](https://consensus.app/papers/details/3c9fe4db5dc7502fa52bdc62f829a2d5/). Soil Tillage Res.
- IPCC, 2019. 2019 Refinement to the 2006 IPCC Guidelines for National Greenhouse Gas Inventories, Vol. 4, Ch. 5 and 11. IPCC, Switzerland.
- Jiang, Z.-W., et al., 2021. [Controlled release urea improves rice production and reduces environmental pollution: a research based on meta-analysis and machine learning](https://consensus.app/papers/details/97387d895e045889af3266bf4f318059/). Environ. Sci. Pollut. Res.
- Lipin, A., et al., 2023. [Modelling nutrient release from controlled release fertilisers](https://consensus.app/papers/details/32cc2243b9a951bbb1b35c66f62d0fd4/). Biosyst. Eng.
- Liu, Y.-Z., et al., 2024. [Localized nitrogen management strategies can halve fertilizer use in Chinese staple crop production](https://consensus.app/papers/details/230e4636223a53d3b02e60251575b41e/). Nat. Food.
- Lu, et al., 2023. **\[TO FILL: full reference for crop-specific CRU effects on NH₃ and N₂O\]**
- Ransom, C.J., et al., 2020. [Nitrogen release rates from slow- and controlled-release fertilizers influenced by placement and temperature](https://consensus.app/papers/details/bf1fcd9b709952e09000b97aeb0ca038/). PLoS ONE.
- Ren, W., et al., 2024. [Evaluating nitrogen dynamic and utilization under controlled-release fertilizer application for sunflowers in an arid region: experimental and modeling approach](https://consensus.app/papers/details/c4b529f3d2c25955919b48fe2020216f/). J. Environ. Manage.
- Shaviv, A., Raban, S., Zaidel, E., 2003. [Modeling controlled nutrient release from polymer coated fertilizers: diffusion release from single granules](https://consensus.app/papers/details/b83f3e5677ad59a8a08500418b797ac2/). Environ. Sci. Technol.
- Verburg, et al., 2020. **\[TO FILL: full reference for soil-moisture effect on PCU release\]**
- Zhang, X., Davidson, E.A., Mauzerall, D.L., Searchinger, T.D., Dumas, P., Shen, Y., 2015. [Managing nitrogen for sustainable development](https://consensus.app/papers/details/831eec19175656b99b0783e37641af27/). Nature.
- Zhang, Y.-J., et al., 2024. [Substituting readily available nitrogen fertilizer with controlled-release nitrogen fertilizer improves crop yield and nitrogen uptake while mitigating environmental risks: a global meta-analysis](https://consensus.app/papers/details/9bcccfc4cdd15ff9af3321e722b8c18b/). Field Crops Res.
- Zhang, et al., 2022. [Can controlled-release urea replace the split application of normal urea in China? A meta-analysis](https://consensus.app/papers/details/36a234488b7150d2a24005c9ade502ed/). Field Crops Res. **\[VERIFY: authors\]**
- Zhao, Z., et al., 2017. [Modelling sugarcane nitrogen uptake patterns to inform design of controlled release fertiliser for synchrony of N supply and demand](https://consensus.app/papers/details/0b347c832efa5b36afeae2dff15c3971/). Field Crops Res.
