# Figure legend — `model-status-graphviz.png`

**Architecture of the platelet Ca²⁺ whole-cell model (v0.4.1).** Schematic
of every process in the calcium pathway, laid out top-to-bottom by
compartment. Agonist input drives three Gq-coupled receptors (P2Y1, PAR1,
PAR4) plus the ATP-gated ionotropic P2X1; the Gαq exchange/GTPase cycle
activates PLCβ, which hydrolyses PIP₂ to IP₃ and DAG; IP₃ opens IP3R
channels in the dense-tubular-system (DTS) membrane, releasing stored Ca²⁺
into the cytosol. SERCA returns Ca²⁺ to the DTS; PMCA and NCX extrude to the
extracellular space; SOCE (STIM1/Orai1) refills the cytosol; mitochondria
buffer transiently via MCU uptake and NCLX release. Cytosolic buffering
(calmodulin, gelsolin) and luminal buffering (calreticulin, HSP90B1, BiP,
CREC) are drawn explicitly. No hand-fitted forcing remains in the calcium
path — every flux is mechanistic.

This is a **schematic, not a result** — it orients the reader before the
data figures; it carries no simulation output.

*Source:* rendered from the Graphviz `dot` block in
`reports/design/model-status-graphviz.qmd` (Quarto). 

*Provenance note:* to be replaced with a hand-drawn BioRender version for the
final manuscript (placeholder schematic for the draft).
