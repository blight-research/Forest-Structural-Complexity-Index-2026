===============================================================================
  FOREST STRUCTURAL COMPLEXITY (FSCI) TOOLSET
  Lidar-based mapping of field-derived forest structural complexity
===============================================================================

OVERVIEW
-------------------------------------------------------------------------------
This toolset lets forest managers and researchers quantify and map forest
structural complexity using airborne lidar. It is based on the framework
developed for the Priest River Experimental Forest (PREF) study in northern
Idaho mixed-conifer forest (Light et al.).

The toolset produces a Forest Structural Complexity Index (FSCI): a continuous,
bounded 0-100 measure built from field-measured structural variables and
predicted across the landscape from lidar.

There are three scripts. Most users need only one or two of them.


WHICH SCRIPT DO I NEED?
-------------------------------------------------------------------------------
1) 01_identify_field_variables.R
   Identifies which field-measured structural variables best distinguish
   canopy-complexity groups. Use this if you are building your own index and
   want to determine which of your field variables to include.
   INPUT : one CSV (plots x lidar complexity metrics x candidate field vars)
   OUTPUT: a table of candidate field variables, flagging those selected.

2) 02A_build_FSCI_model.R
   The full framework. Builds an FSCI index from your chosen field variables,
   then builds and validates a model predicting that index from your lidar
   metrics. Use this if you have BOTH field and lidar data and want a model
   tailored to your own forest.
   INPUT : one CSV (plots x field vars x lidar metrics)
   OUTPUT: a candidate-model table + a saved model file (FSCI_model.rds)

3) 02B_apply_FSCI_model.R
   Applies an existing FSCI model to lidar to produce a wall-to-wall FSCI map.
   Use this if you want to apply the published PREF model (or a model you built
   with 02A) directly to your lidar, without rebuilding anything.
   INPUT : lidar tiles + a DTM + a boundary shapefile + a model (.rds)
   OUTPUT: FSCI map, uncertainty, and a confidence layer (GeoTIFFs)


TYPICAL WORKFLOWS
-------------------------------------------------------------------------------
A) "Just map my forest with the published model"
   -> Run 02B with the provided FSCI_model.rds. (Read the applicability note.)

B) "Build my own index and model, then map"
   -> Run 01  (choose field variables)
   -> Run 02A (build + save your model)
   -> Run 02B (map, using the model you built)


APPLICABILITY (PLEASE READ)
-------------------------------------------------------------------------------
The published model was trained at PREF (northern Idaho mixed-conifer forests). It is
most reliable in structurally and ecologically comparable forests. Script 02B
does NOT hide uncertain predictions; it outputs a confidence layer (Mahalanobis
distance from the training data). Pixels with high values fall outside the
range the model was trained on and should be interpreted with caution. If your
forest differs substantially from PREF, prefer building your own model (02A).


REQUIREMENTS
-------------------------------------------------------------------------------
- R (version 4.0 or newer recommended).
- Each script installs the R packages it needs the first time it runs.
- Script 02B requires lidar tiles, a digital terrain model (DTM), and a
  study-area boundary shapefile, all in the same coordinate system.
- You edit only the clearly marked CONFIG block near the top of each script,
  then run the script from top to bottom.


REPRODUCIBILITY NOTE
-------------------------------------------------------------------------------
Index endpoints: The FSCI index depends on reference "endpoint" values (the
min/max that define 0 and 100). To reproduce the published PREF index exactly,
use the PREF reference endpoints (documented in 02A). To build a relative index
for your own forest, let the script derive endpoints from your own data (the
default).

Map reproduction (Script 2B): Script 2B was validated against the prediction
raster from the original study. The model coefficients are applied exactly, and
two of the three predictor metrics (pz_above_30 and lad_cv) reproduce the
original values to within floating-point precision. A small, consistent
difference remains in the third predictor (zentropy), which propagates to a
median difference of approximately 1.4 FSCI units on the 0-100 scale, with
nearly all cells within about +/- 2.3 units. This residual is most consistent
with a difference in the zentropy computation across versions of the
lidar-processing packages, rather than with the model or the analysis pipeline.
Users reproducing the published map should expect agreement within this small
tolerance rather than bit-for-bit identity. For maximum consistency, record and
match lidar-processing package versions (notably lidR and lidRmetrics).


HOW THESE SCRIPTS WERE PREPARED
-------------------------------------------------------------------------------
These production scripts were adapted from the original research code by the
author with the assistance of Claude (Anthropic). The transformation was guided
by explicit instructions to preserve the original analyses exactly, so that the
cleaned scripts reproduce the same results given the same input data. The author
validated the cleaned scripts against the original outputs.

NOTE- DISCLAIMER
-------------------------------------------------------------------------------

This repository provides an implementation of the structural index described in 
[Light et al]. The code is provided to facilitate application of the index to research and 
forest-management datasets. Users are responsible for determining whether the 
input data, assumptions, and resulting index values are appropriate for their 
particular forest and intended application. The authors recommend independently 
verifying results before using them to inform management decisions.

- the code implements the structural index described in Light et al paper;
- users are responsible for checking that their input data are appropriate;
- users should verify that the code is functioning correctly with their data;
- the authors make no guarantee that the index is appropriate for a particular forest, management objective, or decision;
- the code is provided “as is” and without warranty.


CITATION
-------------------------------------------------------------------------------
If you use this toolset, please cite:
   [Author(s), Year. Title. Journal / repository / DOI.]     <-- In progress...

Contact: Brandon Light aterdennis@gmail.com
License: 
===============================================================================
