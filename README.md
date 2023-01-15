# Flow cytometry analysis on R/bioconductor

This course is a 4 day x 4 hour course instructed by Alice Yue, hosted on Physalia: https://www.physalia-courses.org/courses-workshops/flow-cytometry/

## Course outline

- **Course**: Flow cytometry analysis with R/Bioconductor
- **Date and time**: 2023-01-16 - 2022-01-19 @ 16:00-18:00 CEST (7:00-11:00 PST)
- **Instructor**: Alice Yue (Metafora-biosystems, Paris France, Vancouver Canada)

**Level**: Beginner-Intermediate

**Overview**: Flow cytometry is a gold standard in the analysis of immune cells for clinical diagnosis and research. This course introduces what flow cytometry is and why we use it to analyze cell population composition in biological samples. We will learn about best practices for how to analyze flow cytometry data using R/Bioconductor. Said practices include: preprocessing of the data (compensation, transformation, and quality control), multi-dimensional cell population identification via clustering, and visualizing the results in 2D. These tools can be applied to all types of cytometry including flow, mass, and spectral.

**Target audience**: This course is created for anyone interested in analyzing biological samples with single-cell flow cytometry. Background in flow cytometry and R/Bioconductor is not necessary as we will go over a short introduction to them. However, experience with programming will help greatly.

**Learning outcomes**:
- Understand the flow cytometry machinery and its analytical purpose.
- Be able to set up the infrastructure for and write basic data analytic scripts in R.
- Describe and execute each step in the flow cytometry data analytics pipeline in R/Bioconductor.
- Be comfortable with interpreting and eliciting conclusions from the results of the flow cytometry data analytics pipeline.

**Schedule**:

Support over Slack will be available throughout the course: http://physaliaflowcytometry.slack.com/

| Date       | Session (link to slides) | Script(s) | Practice activities |
|------------|--------------------------|-----------|---------------------|
| 2023-01-16 | 01: [Background on flow cytometry](https://docs.google.com/presentation/d/1O1-l9bhTNjxBxotQkL2kOLpn5keE8tBXuIMzRn0bMw8/edit?usp=sharing) | | |
| 2023-01-16 | 02: [Beginner's guide to R](https://docs.google.com/presentation/d/1PMrVL7BRuhdmD3DsEcP2uYgE3KZIGxtoLDoIkv7s9Bs/edit?usp=sharing) | [01_intro_to_R.R](01_intro_to_R.R), [02_packages](02_packages.R) | Search for "# TRY" in script(s); install the packages in [02_packages](02_packages.R). |
| 2023-01-17 | 03: [Preprocessing flow cytometry data](https://docs.google.com/presentation/d/1HC29MJrkoxpMI59Ezz6yABd1xVn_ZyOjv_SpesqVRKo/edit?usp=sharing) | [03_preprocess.R](03_preprocess.R) | Review the script. |
| 2023-01-18 | 04: [Cell population identification using 2D gating](https://docs.google.com/presentation/d/1dANaxK4I1pclu0IcbQ0UYHcVlsWOk_ni8Otkx_m61gw/edit?usp=sharing) | [04_gating.R](04_gating.R) (Solution: [04_gating_full.R](04_gating_full.R)), [05_gating_stats.R](05_gating_stats.R) | Complete the first script. Review the second script. |
| 2023-01-19 | 05: [Cell population identification using clustering](https://docs.google.com/presentation/d/11bTFEjaQcLnSFef3A8xxEl1a37QcIRJCNrEmvscBUVk/edit?usp=sharing) | [06_gating_clustering](06_gating_clustering.R), [07_clustering.R](07_clustering.R) | Search for "# TRY" in script(s) :) |


**Schedule description**:

- 2023-01-16
    - 01. [Background on flow cytometry](https://docs.google.com/presentation/d/1O1-l9bhTNjxBxotQkL2kOLpn5keE8tBXuIMzRn0bMw8/edit?usp=sharing): This session focuses on the purpose of flow cytometry, a brief overview of its machinery, and its data analysis pipeline.
    - 02. [Beginner's guide to R](https://docs.google.com/presentation/d/1PMrVL7BRuhdmD3DsEcP2uYgE3KZIGxtoLDoIkv7s9Bs/edit?usp=sharing): Users will learn the basics of R and we will set up the computational infrastructure required for the next sessions. For those who do not have R, for this session, please use https://posit.cloud/.
- 2023-01-17
    - 03. [Preprocessing flow cytometry data](https://docs.google.com/presentation/d/1HC29MJrkoxpMI59Ezz6yABd1xVn_ZyOjv_SpesqVRKo/edit?usp=sharing): We will learn about the theory behind preprocessing of flow cytometry samples.
- 2023-01-18
    - 04. [Cell population identification using 2D gating](https://docs.google.com/presentation/d/1dANaxK4I1pclu0IcbQ0UYHcVlsWOk_ni8Otkx_m61gw/edit?usp=sharing): Participants will learn about how flow cytometrists traditionally identify cell populations in preprocessed samples using 2D gating. We will go over an example of an existing manual 2D gating from flowJo and how it is used to compare cell population abundances across samples.
- 2023-01-19
    - 05. [Cell population identification using clustering](https://docs.google.com/presentation/d/11bTFEjaQcLnSFef3A8xxEl1a37QcIRJCNrEmvscBUVk/edit?usp=sharing): We will analyze the same samples using clustering in R and learn how to interpret the results.