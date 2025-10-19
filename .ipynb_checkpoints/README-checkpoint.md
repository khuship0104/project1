# Project 1: Computational Analysis of a Facebook Page-Page Network

**Course:** CAP 6318 - Computational Analysis of Social Complexity  
**Team Members:** Boopathiraj Balasubramani,Kushi Patel



## Project Description

This project analyzes the Facebook Large Page-Page Network from the SNAP dataset to understand its structural properties. Using Julia and the Graphs.jl package, we investigate three core concepts: **homophily** (the tendency for pages to "like" similar pages), the presence of **structural bridges** (pages connecting diverse communities), and **structural balance** (the patterns of triadic closure in the network).



## How to Run the Analysis

This project is built as a reproducible Julia environment.

1.  **Prerequisites:**
    * Julia (v1.12.0 or compatible)
    * A Jupyter environment (like VS Code with the Julia extension, Jupyter Lab, or Jupyter Notebook).

2.  **Setup:**
    * Ensure the `data/` folder, `src/` folder, `Project.toml`, and `Manifest.toml` are all in the same directory as the notebook.
    * Open a Julia REPL (terminal) in this project directory.
    * Run the following commands to install all necessary packages:

    ```julia
    using Pkg
    Pkg.activate(".")
    Pkg.instantiate()
    ```

3.  **Execution:**
    * Open `project1.ipynb` in your Jupyter environment.
    * Run the cells sequentially from top to bottom to replicate the analysis and generate all results and visualizations.



## Key Findings

1.  **Strong Homophily:** The network exhibits strong homophily. Pages are **3.34 times more likely** to connect with other pages from the same category than expected by chance (0.885 observed homophily vs. 0.265 chance). This effect is strongest for "TV Shows" (32.99x) and "Politicians" (11.64x).

2.  **Bridging Nodes:** Structural bridges, which connect different communities, are predominantly **Companies** and **Government Organizations**. These page types play a key role in connecting otherwise disparate parts of the network, such as linking politician, company, and TV show clusters.

3.  **Structural Balance:** The network shows clear signs of structural balance. We found that **99.9% of all triangles are balanced** (all positive links). Furthermore, the "friend-of-a-friend" closure rate is **26.18%**, which is **1.13 times higher** than the random baseline (23.26%), indicating that if two pages like a common third page, they are more likely to like each other as well.