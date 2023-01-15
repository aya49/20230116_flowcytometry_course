# Flow cytometry analysis on R/bioconductor

This course is a 4 day x 4 hour course taught by Alice Yue, hosted on Physalia: https://www.physalia-courses.org/courses-workshops/flow-cytometry/

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

**Program**:
- Support over Slack will be available throughout the course: <link>
- 2023-01-16
    - [Session 1 Background on flow cytometry](https://docs.google.com/presentation/d/1O1-l9bhTNjxBxotQkL2kOLpn5keE8tBXuIMzRn0bMw8/edit?usp=sharing): This session will focus on the purpose of flow cytometry, a brief overview of its machinery, and what the data output looks like. We will also go over the steps in the flow cytometry data analysis pipeline.
    - [Session 2 Beginner's guide to R](https://docs.google.com/presentation/d/1PMrVL7BRuhdmD3DsEcP2uYgE3KZIGxtoLDoIkv7s9Bs/edit?usp=sharing) [scripts: [01_intro_to_R.R](./01_intro_to_R.R), [02_packages.R](./02_packages.R)]: Users will learn about the basics of R and we will set up the computational infrastructure for R on https://rstudio.cloud/.
- 2023-01-17
    - Session 3 Preprocessing flow cytometry data: We will go through a theoretical overview of the steps in preprocessing flow cytometry samples. A practical hands-on walkthrough in R will follow.
- 2023-01-18
    - Session 4 Cell population identification using 2D gating: Participants will learn about how flow cytometrists traditionally identify cell populations in preprocessed samples using 2D gating. We will go over an example of an existing manual 2D gating from flowJo and how it is used to compare cell population abundances across samples. We will also go over how to replicate manual gating in R with flowDensity.
- 2023-01-19
    - Session 5 Cell population identification using clustering: We will analyze the same samples using clustering in R and learn how to interpret the results.
