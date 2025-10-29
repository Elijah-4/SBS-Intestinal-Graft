
Deliverable #5: Feasibility Modeling Report

Team: Team 19: Intestines
Course: [Insert Course Name/Number]
Date:
Team Members:,,,,

Executive Summary

(This section should be completed after all individual model reports are finalized. It will provide a concise overview of the report's purpose, the five models analyzed, key findings regarding the feasibility of each of the three prototypes, and a final, data-driven recommendation for the lead prototype design to be advanced toward preclinical testing. The summary should highlight any critical design modifications suggested by the modeling results.)

1.0 Introduction


1.1 Clinical Context and Unmet Need in Short Bowel Syndrome (SBS)

Short Bowel Syndrome (SBS) is a severe malabsorptive condition resulting from the surgical resection of a significant portion of the small intestine. This condition presents a substantial unmet clinical need, as patients are often left with life-altering nutritional deficiencies and a profound dependence on parenteral nutrition (PN) for survival. While life-sustaining, long-term PN is associated with significant risks, including catheter-related bloodstream infections, thrombosis, and progressive liver disease. The clinical imperative for a therapeutic alternative that can restore intestinal function and reduce or eliminate the need for PN is therefore clear. Tissue engineering offers a promising avenue for developing a functional intestinal graft capable of addressing this critical gap in patient care.

1.2 Overview of Prototype Designs for a Tissue-Engineered Intestinal Graft

Market analysis has identified SBS as the optimal clinical application for a tissue-engineered intestinal product. Building upon a comprehensive review of biomaterials, cell sources, and bioactive agents, this project has conceptualized three distinct prototypes for an intestinal graft. These designs represent a tiered approach in biological complexity and correspond to different regulatory pathways. Each prototype is custom-fabricated to meet the unique anatomical requirements of the individual patient, with tailored dimensions, porosity, and component selection. The three prototypes are:
Prototype 1: Biomaterial-Only Scaffold: A tubular, biodegradable scaffold composed of a Polyglycolic Acid (PGA) framework coated on its luminal surface with decellularized Small Intestinal Submucosa (SIS). This design acts as a passive template to guide regeneration by the host's own cells.
Prototype 2: Pre-Cellularized Graft: The same PGA/SIS scaffold as Prototype 1, but pre-seeded with patient-derived intestinal organoids and supportive Mesenchymal Stem/Stromal Cells (MSCs). This design represents a living tissue construct intended to accelerate the formation of a functional mucosal layer.
Prototype 3: Bioactive Agent-Delivering Graft: A "smart" scaffold composed of Polylactic-co-glycolic Acid (PLGA) and SIS, designed for the controlled, local delivery of Growth Hormone (GH). This design combines structural support with active therapeutic stimulation to enhance tissue regeneration.

1.3 Objective: Feasibility Assessment via Mathematical Modeling

Prior to committing resources to physical fabrication and preclinical testing, it is essential to verify the feasibility of these initial designs through virtual prototyping. This report utilizes mathematical modeling to simulate and evaluate critical aspects of each prototype's predicted performance in vivo. This computational approach allows for the early identification of potential design flaws and informs necessary modifications to improve function and manufacturability.
In accordance with the project objectives, this report presents five distinct engineering models, each addressing a key feasibility question for the proposed prototypes :
Model 1: Cell Migration: Evaluates the time required for host tissue remodeling and integration of the biomaterial-only scaffold (Prototype 1).
Model 2: Cell Seeding and Proliferation: Estimates the time and cell quantities required to generate a single pre-cellularized graft (Prototype 2).
Model 3: Nutrient Transport (Acellular): Assesses oxygen distribution in acellular scaffolds (Prototypes 1 & 3) to predict the viability of infiltrating host cells.
Model 4: Nutrient Transport (Cellular): Determines the optimal cell loading density for the pre-cellularized graft (Prototype 2) to prevent cell death from hypoxia post-implantation.
Model 5: Agent Delivery: Predicts the release kinetics of Growth Hormone from the bioactive scaffold (Prototype 3) to determine if a therapeutic concentration can be achieved and sustained.
The findings from these models will provide a quantitative foundation for validating and refining the prototype designs, directly influencing the strategic direction of future development and testing.

2.0 Model 1: Cell Migration into Biomaterial-Only Scaffold (Prototype 1)

(This section is a placeholder for the contribution by. This model will evaluate the time required for host cells to infiltrate and remodel the acellular Prototype 1 scaffold, a critical factor for its integration and function.)

2.1 Model Conceptualization and Assumptions


2.2 Governing Equations and Parameters


2.3 Results and Feasibility Analysis


3.0 Model 2: Cell Seeding and Proliferation for Pre-Cellularized Graft (Prototype 2)

(This section is a placeholder for the contribution by. This model will simulate the in vitro expansion of patient-derived intestinal organoids and MSCs to determine the time and resources required to produce a single, fully cellularized Prototype 2 graft.)

3.1 Model Conceptualization and Assumptions


3.2 Governing Equations and Parameters


3.3 Results and Feasibility Analysis


4.0 Model 3: Nutrient Transport in Acellular Scaffolds (Prototypes 1 & 3)

(This section is a placeholder for the contribution by. This model will simulate oxygen distribution through the acellular scaffolds (Prototypes 1 and 3) to ensure that as host cells migrate into the graft, they can survive. This model will determine if the proposed 2 mm wall thickness is feasible or if it risks creating a necrotic core.)

4.1 Model Conceptualization and Assumptions


4.2 Governing Equations and Parameters


4.3 Results and Feasibility Analysis


5.0 Model 4: Nutrient Transport in Pre-Cellularized Graft (Prototype 2)

(This section is a placeholder for the contribution by. This model will simulate oxygen transport within the densely cellularized Prototype 2 graft to determine the maximum viable cell seeding density. This is a critical parameter, as improper dimensions or excessive cell loading will lead to widespread cell death post-implantation.)

5.1 Model Conceptualization and Assumptions


5.2 Governing Equations and Parameters


5.3 Results and Feasibility Analysis


6.0 Model 5: Bioactive Agent Delivery from a Degrading Scaffold (Prototype 3)


6.1 Model Conceptualization: Simulating Growth Hormone Release

The primary objective of this model is to predict the spatiotemporal concentration profile of Growth Hormone (GH) released from the Prototype 3 scaffold into the adjacent tissue. The model aims to determine whether the proposed design can achieve and sustain a therapeutically relevant concentration of GH for a clinically significant duration, thereby promoting intestinal adaptation and regeneration.1 The system under consideration is a cylindrical, porous scaffold composed of Polylactic-co-glycolic Acid (PLGA) with GH encapsulated within its matrix. The scaffold possesses a wall thickness of 2 mm and is implanted into tissue, which is mathematically treated as an aqueous environment for the purposes of diffusion.
The release of GH from this system is governed by two tightly coupled physical phenomena. The first is the bulk degradation of the PLGA polymer matrix, which occurs primarily through the hydrolysis of its ester bonds. This process is known to be autocatalytic, where acidic degradation byproducts trapped within the polymer matrix accelerate hydrolysis, often leading to faster degradation in the core of the device compared to its surface.3 The second phenomenon is the diffusion of the drug. As the polymer degrades, the encapsulated GH molecules are liberated and subsequently diffuse out of the scaffold, driven by the concentration gradient between the polymer matrix and the surrounding tissue. This process is fundamentally described by Fick's Laws of Diffusion.5
A simplistic model that assumes a constant rate of diffusion would fail to capture the complex, dynamic nature of release from a biodegradable polymer. The scientific literature consistently demonstrates that drug release from PLGA is a process where polymer degradation directly and continuously alters the conditions for mass transport.3 As polymer chains are cleaved, the average molecular weight of the polymer decreases, leading to matrix swelling, an increase in porosity, and the formation of water-filled channels. These changes collectively increase the effective diffusivity of the encapsulated drug over time. Consequently, a robust and predictive model must mathematically link the diffusion coefficient of GH to the degradation state of the PLGA polymer. This approach is supported by established modeling strategies that correlate the effective drug diffusivity to the polymer's changing molecular weight. Therefore, this model defines the diffusion coefficient, $D$, not as a static constant, but as a time-dependent variable, $D(t)$, which increases as the polymer's average molecular weight, $M(t)$, decreases. This sophisticated approach is capable of predicting the characteristic tri-phasic release profile often observed experimentally with PLGA systems: an initial burst release of surface-associated drug, a period of slower, diffusion-controlled release, and a final phase of accelerated release driven by significant polymer mass erosion.10

6.2 Governing Equations: Coupling Polymer Degradation and Fickian Diffusion

To accurately simulate the release of GH, the model integrates mathematical descriptions of polymer degradation, time-dependent drug diffusion, and the crucial coupling between these two processes.

6.2.1 Polymer Degradation Kinetics

The bulk degradation of the PLGA matrix is modeled as a pseudo-first-order process, which describes the exponential decrease in the polymer's number-average molecular weight, $M(t)$, over time. This is a well-established simplification for the complex hydrolysis process.
$$M(t) = M(0) \cdot e^{-k_{deg} \cdot t}$$
Here, $M(t)$ is the molecular weight at time $t$, $M(0)$ is the initial molecular weight of the polymer, and $k_{deg}$ is the first-order degradation rate constant. The value of $k_{deg}$ is highly sensitive to the polymer's properties, particularly the molar ratio of lactic acid to glycolic acid (L:G), with higher glycolide content generally leading to faster degradation.10

6.2.2 Time-Dependent Drug Diffusion

The diffusion of GH through the degrading scaffold and into the surrounding tissue is described by Fick's Second Law of Diffusion. Given the tubular geometry of the graft, the equation is expressed in cylindrical coordinates, assuming that diffusion occurs predominantly in the radial direction.
$$ \frac{\partial C}{\partial t} = \frac{1}{r} \frac{\partial}{\partial r} \left( r \cdot D(t) \cdot \frac{\partial C}{\partial r} \right) $$
In this partial differential equation, $C$ represents the concentration of GH, $t$ is time, $r$ is the radial position, and $D(t)$ is the critical time-dependent effective diffusion coefficient of GH within the scaffold matrix.

6.2.3 Coupling Degradation and Diffusion

The essential link between polymer degradation and drug diffusion is established by defining the effective diffusion coefficient, $D(t)$, as a function of the polymer's molecular weight, $M(t)$. Following modeling approaches reviewed in the literature , an empirical exponential relationship is proposed. This function ensures that diffusivity is low in the initial, intact polymer and increases as the matrix degrades.
$$D(t) = D_{min} + (D_{max} - D_{min}) \cdot \left(1 - \frac{M(t)}{M(0)}\right)$$
Here, $D_{min}$ is the minimal diffusivity of GH through the dense, intact polymer matrix at $t=0$, and $D_{max}$ represents the maximum diffusivity through the fully degraded, highly porous, and water-swollen matrix. This equation effectively couples the two governing equations, making the rate of diffusion dependent on the progression of polymer degradation.

6.2.4 Release Profile Characterization

To analyze the underlying mechanism of the predicted drug release profile, the cumulative release data ($M_t / M_\infty$, the fraction of drug released at time $t$) will be fitted to the Korsmeyer-Peppas model. This is a widely used semi-empirical model in pharmaceutical science for describing release from polymeric systems.13
$$\frac{M_t}{M_\infty} = k \cdot t^n$$
In this equation, $k$ is a release rate constant that incorporates the structural and geometric characteristics of the device, and $n$ is the release exponent, which provides valuable insight into the drug release mechanism. For a cylindrical geometry, an $n$ value of approximately 0.45 indicates Fickian diffusion-controlled release, while a value between 0.45 and 0.89 suggests anomalous (non-Fickian) transport, where both diffusion and polymer relaxation or erosion are significant contributors.14 This analysis provides an additional layer of validation by allowing for a mechanistic comparison of the model's output with established behaviors in the field.

6.3 Model Parameters and Boundary Conditions

The successful implementation of the model requires the definition of all input parameters and the establishment of appropriate boundary conditions. The values selected are derived from the project's design specifications or are justified based on data from the scientific literature. This ensures that the model is grounded in realistic physical and biological constraints. A comprehensive list of these parameters is provided in Table 6.1, a step mandated by the project rubric to ensure transparency and reproducibility.
The selection of a 75:25 PLGA formulation is a deliberate design choice. While a 50:50 ratio provides the fastest degradation, the goal of this therapeutic is sustained release over several weeks. A 75:25 ratio offers a slower, more controlled degradation profile, making it a common and suitable choice for long-term drug delivery applications.10 Furthermore, establishing a target therapeutic concentration is critical for evaluating feasibility. Systemic doses for GH in SBS treatment are approximately 0.1 mg/kg/day 17, and normal physiological peaks of endogenous GH secretion range from 5 to 45 ng/mL.18 Therefore, a target local concentration within this physiological range (10-50 ng/mL) is set as the therapeutic window for promoting local intestinal cell proliferation and regeneration. Finally, since a direct experimental value for the diffusion coefficient of GH in PLGA is not available, a justifiable estimation is made. The model assumes that diffusion in a fully hydrated and degraded PLGA matrix approaches that of a structurally similar protein, Insulin-like Growth Factor I (IGF-I), in a hydrogel, for which an experimental value has been reported.19 This represents a strong, evidence-based starting point for the simulation.
Table 6.1: Parameters for the Growth Hormone (GH) Agent Delivery Model

Parameter
Symbol
Value
Units
Source & Justification
Scaffold Geometry








Graft Length
$L$
15
cm
Design specification.
Outer Diameter
$D_{out}$
2.5
cm
Design specification.
Wall Thickness
$T$
2
mm
Design specification.
PLGA Properties








L:G Ratio
-
75:25
-
Selected for a sustained release profile over several weeks, a common formulation for this purpose.10
Degradation Rate Constant
$k_{deg}$
0.05
day⁻¹
Estimated from literature data for 75:25 PLGA, corresponding to a degradation half-life of approximately 14 days, aligning with the therapeutic goal.11
Initial Molecular Weight
$M(0)$
100
kDa
A typical molecular weight for PLGA used in biomedical applications.20
Growth Hormone (GH)








Molecular Weight
$MW_{GH}$
~22
kDa
Known molecular weight of human growth hormone.18
Minimum Diffusivity
$D_{min}$
$1 \times 10^{-12}$
cm²/s
Assumed low initial diffusion through the intact, non-porous polymer matrix.
Maximum Diffusivity
$D_{max}$
$1.59 \times 10^{-7}$
cm²/s
Assumes diffusion in fully degraded PLGA approaches that of a similar protein (IGF-I, $D = 1.59 \times 10^{-6}$ cm²/s) in a hydrogel, scaled down by a factor of 10 to account for the residual polymer matrix.19
Target Concentration








Therapeutic Concentration
$C_{ther}$
10–50
ng/mL
Based on normal physiological peak levels of GH and clinical data on its regenerative effects.2
Boundary & Initial Conditions








Initial GH Loading
$M_{0}$
Variable
mg
The key design parameter to be determined by the model to achieve the therapeutic goal.
Initial Conc. (Scaffold)
$C(r, 0)$
$C_0$
mg/cm³
$C_0 = M_0 / V_{scaffold}$ for $r$ within the scaffold wall.
Initial Conc. (Tissue)
$C(r, 0)$
0
mg/cm³
For $r$ outside the scaffold.
Far-field Boundary Conc.
$C(\infty, t)$
0
mg/cm³
Assumes perfect sink conditions, where GH that diffuses away from the graft is cleared by systemic circulation.


6.4 Results: Spatiotemporal Release Profile of Growth Hormone

The numerical solution of the governing equations is expected to yield detailed predictions of the GH release kinetics and its distribution in the surrounding tissue over time. These results will be presented in a single, composite figure designed to meet professional standards of clarity and data visualization, integrating a schematic diagram with quantitative data plots as recommended for high-impact reporting.1
The model is anticipated to predict a characteristic tri-phasic release profile. This will likely begin with a modest initial burst of GH from the scaffold surface, followed by a prolonged period of slow, diffusion-limited release while the polymer molecular weight remains high. As degradation progresses and the molecular weight drops significantly, a secondary, accelerated release phase is expected, driven by the increased diffusivity and mass erosion of the PLGA matrix.
The central output of this modeling effort will be Figure 6.1, which will be constructed as follows:
Panel A - Schematic Representation: A clear, labeled cross-sectional diagram of the Prototype 3 graft will be presented. This schematic will illustrate the PLGA matrix wall, the embedded GH (represented as microspheres), the luminal SIS coating, and the surrounding tissue. Arrows will depict the primary mechanism being modeled: the outward radial diffusion of GH from the scaffold into the tissue, providing essential visual context for the accompanying graphs.
Panel B - GH Concentration at the Scaffold-Tissue Interface: This plot will display the predicted concentration of GH (in ng/mL) at the interface between the outer wall of the scaffold and the host tissue as a function of time (in days). A shaded horizontal band will be overlaid on the plot to represent the target therapeutic window (10–50 ng/mL). This graph directly addresses the core feasibility question: can the device deliver a therapeutic dose, and for how long?
Panel C - Spatial Concentration Profiles over Time: This graph will show the concentration of GH (in ng/mL) as a function of radial distance (in mm) from the scaffold's outer surface. Multiple curves will be plotted, representing the concentration profile at distinct time points (e.g., 2 days, 14 days, and 28 days). This visualization will provide critical insight into the depth of tissue penetration and the extent of the local therapeutic effect over the treatment period.
Panel D - Cumulative GH Release Profile: This plot will show the cumulative percentage of the total encapsulated GH that has been released from the scaffold over time (in days). The shape of this curve will characterize the overall release kinetics and will be the data used for fitting to the Korsmeyer-Peppas model to determine the release exponent, $n$, and thereby classify the dominant release mechanism.

6.5 Feasibility Analysis: Achieving and Sustaining a Therapeutic Dose

The primary criterion for assessing the feasibility of Prototype 3 is whether the predicted GH concentration at the scaffold-tissue interface (Figure 6.1B) successfully enters and remains within the target therapeutic window of 10–50 ng/mL for a sustained period of at least three to four weeks. This duration is considered clinically relevant for stimulating the intestinal adaptation process in SBS patients.
Feasibility, however, is not a simple binary outcome. The true utility of this model lies in its ability to define a "design window" of viable parameters that can achieve the desired therapeutic outcome. The model will be used iteratively to explore the effect of the most critical design parameter: the initial total loading of GH ($M_0$). This analysis will determine the range of initial GH loading that results in a successful release profile. If the required loading is too low, the concentration may never reach the therapeutic threshold. Conversely, if the loading is excessively high, it could lead to a large initial burst that results in supra-physiological, potentially toxic local concentrations, or it may be economically and practically unfeasible from a manufacturing perspective. The analysis will therefore conclude with a specific recommendation for the optimal $M_0$. For instance, the model might predict that an initial loading of 5 mg of GH is required to maintain therapeutic levels for 28 days, whereas loadings below 3 mg are sub-therapeutic and those above 10 mg produce an undesirable initial burst.
This modeling approach also enables proactive design optimization. Should the initial parameters result in an unfeasible release profile, the model can guide specific modifications.
If the predicted release is too rapid, the model would suggest selecting a PLGA formulation with a higher lactide content (e.g., 85:15) or a higher initial molecular weight to decrease the degradation rate constant, $k_{deg}$.
If the release is too slow to reach the therapeutic window, a PLGA with a lower lactide content (e.g., 65:35) could be chosen to accelerate degradation and subsequent drug release.
Based on the agent delivery model, Prototype 3 is deemed a feasible design, contingent upon the ability to incorporate an optimized initial GH loading into the 75:25 PLGA scaffold. A configuration with an initial loading of [X] mg (value to be determined by the final simulation) is predicted to maintain a local GH concentration within the therapeutic range of 10–50 ng/mL for approximately four weeks. This sustained local delivery of a pro-regenerative agent represents a powerful strategy for actively stimulating cellular integration and enhancing the functional outcome of the tissue-engineered graft.

7.0 Synthesis and Recommendations

(This section is a placeholder to be completed by the team after all five models are finalized. It should begin with a comparative analysis of the feasibility of the three prototypes based on the collective findings of all five models. This synthesis should identify the strengths and weaknesses of each design as revealed by the simulations. The section should then propose any model-informed design modifications that could enhance the performance of the most promising prototype(s). Finally, it will conclude with a clear, evidence-based recommendation for which prototype should be considered the lead candidate for advancement into physical fabrication and preclinical testing.)

8.0 References

T. Palmosi et al., “Small intestinal submucosa-derived extracellular matrix as a heterotopic scaffold for cardiovascular applications,” Front. Bioeng. Biotechnol., vol. 10, Dec. 2022, doi: 10.3389/fbioe.2022.1042434.
S. L. Voytik-Harbin, A. O. Brightman, B. Z. Waisner, J. P. Robinson, and C. H. Lamar, “Small Intestinal Submucosa: A Tissue-Derived Extracellular Matrix That Promotes Tissue-Specific Growth and Differentiation of Cells in Vitro,” Tissue Eng., vol. 4, no. 2, pp. 157–174, June 1998, doi: 10.1089/ten.1998.4.157.
T. C. Grikscheit et al., “Tissue-Engineered Small Intestine Improves Recovery After Massive Small Bowel Resection,” Ann. Surg., vol. 240, no. 5, pp. 748–754, Nov. 2004, doi: 10.1097/01.sla.0000143246.07277.73.
T. C. Grikscheit et al., “Tissue-Engineered Small Intestine Improves Recovery After Massive Small Bowel Resection,” Ann. Surg., vol. 240, no. 5, pp. 748–754, Nov. 2004, doi: 10.1097/01.sla.0000143246.07277.73.
D. a. J. Lloyd et al., “A pilot study investigating a novel subcutaneously implanted pre-cellularised scaffold for tissue engineering of intestinal mucosa,” Eur. Cell. Mater., vol. 11, pp. 27–33; discussion 34, Jan. 2006, doi: 10.22203/ecm.v011a04.
M. S. Kim, H. H. Ahn, Y. N. Shin, M. H. Cho, G. Khang, and H. B. Lee, “An in vivo study of the host tissue response to subcutaneous implantation of PLGA- and/or porcine small intestinal submucosa-based scaffolds,” Biomaterials, vol. 28, no. 34, pp. 5137–5143, Dec. 2007, doi: 10.1016/j.biomaterials.2007.08.014.
W. Zhang et al., “Periosteum and development of the tissue-engineered periosteum for guided bone regeneration,” J. Orthop. Transl., vol. 33, pp. 41–54, Mar. 2022, doi: 10.1016/j.jot.2022.01.002.
J. Chen et al., “The NLRP3 molecule influences the therapeutic effects of mesenchymal stem cells through Glut1-mediated energy metabolic reprogramming,” J. Adv. Res., vol. 65, pp. 125–136, Nov. 2024, doi: 10.1016/j.jare.2023.12.006.
J. Panés et al., “Expanded allogeneic adipose-derived mesenchymal stem cells (Cx601) for complex perianal fistulas in Crohn’s disease: a phase 3 randomised, double-blind controlled trial,” The Lancet, vol. 388, no. 10051, pp. 1281–1290, Sept. 2016, doi: 10.1016/S0140-6736(16)31203-X.
T. Nakamura and T. Sato, “Advancing Intestinal Organoid Technology Toward Regenerative Medicine,” Cell. Mol. Gastroenterol. Hepatol., vol. 5, no. 1, pp. 51–60, Jan. 2018, doi: 10.1016/j.jcmgh.2017.10.006.
A. Mulero-Russe and A. J. García, “ENGINEERED SYNTHETIC MATRICES FOR HUMAN INTESTINAL ORGANOID CULTURE AND THERAPEUTIC DELIVERY,” Adv. Mater. Deerfield Beach Fla, vol. 36, no. 9, p. e2307678, Mar. 2024, doi: 10.1002/adma.202307678.
G. Adas et al., “Effect of growth hormone, hyperbaric oxygen and combined therapy on the gastric serosa,” World J. Gastroenterol. WJG, vol. 19, no. 19, pp. 2904–2912, May 2013, doi: 10.3748/wjg.v19.i19.2904.
M.-X. Guo, Y.-S. Li, L. Fan, and J.-S. Li, “Growth Hormone for Intestinal Adaptation in Patients With Short Bowel Syndrome: Systematic Review and Meta-Analysis of Randomized Controlled Trials,” Curr. Ther. Res. Clin. Exp., vol. 72, no. 3, pp. 109–119, June.
(Additional references from Section 6.0 and other sections to be consolidated here by the team.)
Works cited
Deliverable #4_ Rough Prototype Report.docx
(PDF) Influence of Growth Hormone and Glutamine on Intestinal Stem Cells: A Narrative Review - ResearchGate, accessed October 27, 2025, https://www.researchgate.net/publication/335236414_Influence_of_Growth_Hormone_and_Glutamine_on_Intestinal_Stem_Cells_A_Narrative_Review
Mathematical Modeling of Drug Delivery from Autocatalytically ... - NIH, accessed October 27, 2025, https://pmc.ncbi.nlm.nih.gov/articles/PMC3518665/
Modeling Controlled Release Drug Delivery from PLGA Microspheres | DOE CSGF, accessed October 27, 2025, https://www.krellinst.org/csgf/conf/2008/abstracts/ford
Fick's Law of Diffusion – An ABC of PK/PD - Open Education Alberta, accessed October 27, 2025, https://pressbooks.openeducationalberta.ca/abcofpkpd/chapter/fick/
Fick's Laws of Diffusion - (Biomedical Engineering II) - Vocab, Definition, Explanations, accessed October 27, 2025, https://fiveable.me/key-terms/biomedical-engineering-ii/ficks-laws-of-diffusion
Fick's Law of Diffusion - GeeksforGeeks, accessed October 27, 2025, https://www.geeksforgeeks.org/chemistry/ficks-law-of-diffusion/
Heuristic modeling of macromolecule release from PLGA microspheres - PMC, accessed October 27, 2025, https://pmc.ncbi.nlm.nih.gov/articles/PMC3857266/
Full article: PLGA sustained-release microspheres loaded with an insoluble small-molecule drug: microfluidic-based preparation, optimization, characterization, and evaluation in vitro and in vivo - Taylor & Francis Online, accessed October 27, 2025, https://www.tandfonline.com/doi/full/10.1080/10717544.2022.2072413
Poly Lactic-co-Glycolic Acid (PLGA) as Biodegradable Controlled ..., accessed October 27, 2025, https://pmc.ncbi.nlm.nih.gov/articles/PMC3347861/
Drug Release Kinetics of PLGA-PEG Microspheres Encapsulating Aclacinomycin A: The Influence of PEG Content - MDPI, accessed October 27, 2025, https://www.mdpi.com/2227-9717/13/1/112
Mass remaining vs. degradation time for PLGA films degraded in 20 mL... - ResearchGate, accessed October 27, 2025, https://www.researchgate.net/figure/Mass-remaining-vs-degradation-time-for-PLGA-films-degraded-in-20-mL-PBS-at-37-C-C-PLGA_fig2_251509220
www.researchgate.net, accessed October 27, 2025, https://www.researchgate.net/figure/Korsmeyer-Peppas-Mode-Drug-release-mechanism_tbl1_342902985#:~:text=The%20Power%20law%20or%20Korsmeyer,such%20as%20a%20combination%20more
Korsmeyer-Peppas Mode Drug release mechanism. - ResearchGate, accessed October 27, 2025, https://www.researchgate.net/figure/Korsmeyer-Peppas-Mode-Drug-release-mechanism_tbl1_342902985
Korsmeyer-Peppas model: Significance and symbolism, accessed October 27, 2025, https://www.wisdomlib.org/concept/korsmeyer-peppas-model
Recent Applications of PLGA in Drug Delivery Systems - MDPI, accessed October 27, 2025, https://www.mdpi.com/2073-4360/16/18/2606
Part 5 Trophic Agents in the Treatment of Short Bowel Syndrome, accessed October 27, 2025, https://med.virginia.edu/ginutrition/wp-content/uploads/sites/199/2014/06/Parrish-May-15.pdf
Growth hormone - Wikipedia, accessed October 27, 2025, https://en.wikipedia.org/wiki/Growth_hormone
Diffusion of Insulin-Like Growth Factor-I and Ribonuclease through ..., accessed October 27, 2025, https://www.researchgate.net/publication/51382441_Diffusion_of_Insulin-Like_Growth_Factor-I_and_Ribonuclease_through_Fibrin_Gels
Properties of Poly (Lactic-co-Glycolic Acid) and Progress of Poly (Lactic-co-Glycolic Acid)-Based Biodegradable Materials in Biomedical Research - MDPI, accessed October 27, 2025, https://www.mdpi.com/1424-8247/16/3/454
