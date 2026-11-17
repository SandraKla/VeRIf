####################################### Texts #####################################################
# Use Unicode escapes for non-ASCII characters so sourcing also works in the C locale.

data_text <- HTML(paste0(
  "Load demo data, upload your own file, or enter values in the editable table below. ",
  "Analysis tabs become available once data are loaded. Then use the sidebar to select the laboratory analyte, category, sex and age range."
))
upload_data_text <- HTML(paste0(
  "<strong>File requirements:</strong> Include column headers, age in years, sex and at least one laboratory analyte (use its name as the column header). Columns may be in any order. ",
  "CSV files must use semicolons (;) between values and commas (,) for decimals.", br(),
  "After upload, check the age, sex and optional category assignments and the female/male values. Choose None (all rows) if no category is needed. ",
  "Uploaded files are applied automatically; Submit is only needed for the editable table."
))
reflim_text <- HTML(paste0(
  "These tab displays the corresponding plot and the outputs of the reflim() function, providing an estimation of new reference intervals or a verification of the selected target values. 
  EL stands for equivalence limit. UM stands for uncertainty margin.
  By clicking Visualization of all plots across every process step, all plots generated throughout the workflow can be displayed."
))
refineR_text<- HTML(paste0(
  "If, during the verification with reflimR and its target values or own target values, a yellow or red bar appears, 
  a follow-up analysis using refineR is recommended. The resulting reference intervals from refineR can be used as new target values
  and re-verified with reflimR. If all indicators turn green, this suggests that the manufacturer\u2019s target values are likely incorrect.
  If one or more indicators remain yellow or red, the data are considered too challenging for indirect methods. This assumption can be further 
  evaluated in the mclust tab using a Gaussian mixture model (mclust)."
))
mclust_text <- HTML(paste0(
  "Gaussian mixture modelling for the verification of reference intervals. ",
  "Select or type a cluster count (1-15), or choose Automatic (BIC) to select the count with the highest BIC."
))
mclust_bic_text <- HTML(paste0(
  "Each point shows the Bayesian Information Criterion (BIC) for a candidate number of clusters, using the selected variance model and distribution. ",
  "Higher BIC values indicate a better balance between fit and model complexity; the highest point is preferred among the fitted candidates. ",
  "Similar BIC values indicate a less clear preference. ",
  "This plot compares 1-15 clusters even when a fixed number is selected for the mixture plot."
))
scatterplot_text <- HTML(paste0(
  "The plots show the relationship between age, sex and the laboratory value."
))
zlog_text <- HTML(paste0(
  "zlog values are calculated from the dataset and the calculated reference intervals. The lower reference limits (LL) and upper reference limits (UL) can transform any result x into a zlog value using the following equation:
  zlog(x) = (log(x)\u2013(log(LL)+ log(UL))/2)*3.92/(log(UL)\u2013log(LL)). Values ranging from \u20131.96 to 1.96 are considered normal, while values below \u20135 and above 5 indicate pathological conditions."
))
rpart_text <- HTML(paste0(
  "Decision trees are used here to divide the data into subgroups based on age and sex. The method identifies groups with similar values, allowing more appropriate calculation of reference intervals. The tree is built using the rpart package and visualized with rpart.plot."
))
about_text <- HTML(paste0(
  "<p>VeRIf is an interactive Shiny web application for the verification and evaluation of reference intervals based on routine laboratory data using the R package reflimR. The web application supports medical laboratories in efficiently, transparently, and data-driven reviewing existing reference intervals.</p>",
  br(), "<table class='table table-condensed' style='width: 100%; max-width: 900px;'>",
  "<tr><td><strong>Version:</strong></td><td>1.0.6</td></tr>",
  "<tr><td><strong>Depends:</strong></td><td>R (&gt;= 4.5.2)</td></tr>",
  "<tr><td><strong>Imports:</strong></td><td>DT, mclust, refineR, reflimR, rhandsontable, readxl, rpart, rpart.plot, shiny, shinycssloaders, shinydashboard</td></tr>",
  "<tr><td><strong>Author:</strong></td><td>Sandra Klawitter</td></tr>",
  "<tr><td><strong>E-Mail:</strong></td><td>sandra.klawitter@trillium.de</td></tr>",
  "<tr><td><strong>Workshop:</strong></td><td><a href='https://github.com/SandraKla/VeRIf-workshop'>https://github.com/SandraKla/VeRIf-workshop</a></td></tr>",
  "<tr><td><strong>BugReports:</strong></td><td><a href='https://github.com/SandraKla/VeRIf/issues'>https://github.com/SandraKla/VeRIf/issues</a></td></tr>",
  "<tr><td><strong>R-Code:</strong></td><td><a href='https://github.com/SandraKla/VeRIf'>https://github.com/SandraKla/VeRIf</a></td></tr>",
  "</table>", br(),
  "<p><strong>Disclaimer:</strong> Only anonymized data may be uploaded to this application. This application is provided &quot;as is&quot; and &quot;as available&quot;, without any express or implied warranties of any kind. No warranty is given regarding the accuracy, completeness, reliability, or timeliness of the results. ",
  "The results are provided for informational and research purposes only and must not be used for diagnosis, treatment, prevention, or any form of clinical or medical decision-making. This application is not a medical device or medical product and does not replace professional medical advice. ",
  "To the fullest extent permitted by law, the author disclaims all liability for any direct, indirect, incidental, consequential, or special damages arising from the use of this application or its results. Use of this application is entirely at your own risk.</p>"
))
