# Gene Expression Analysis

## Overview

`gene-expression-analysis` provides Python tools for analyzing gene expression data. The repository includes methods for calculating Pearson correlation coefficients between genes, visualizing the distribution of these correlations, and identifying key genes based on their correlation with a specific gene, MCTS1.

## Features

- **Pearson Correlation Calculation:** Compute Pearson correlation coefficients for all pairs of genes in your dataset.
- **Correlation Distribution Visualization:** Plot the distribution of correlation coefficients to understand the overall gene correlation behavior.
- **Key Gene Analysis:** Identify genes with the highest positive, negative, and least correlations with MCTS1.
- **Scatter Plots with Regression:** Create scatter plots with linear regression fits to explore relationships between MCTS1 and other genes.

## Getting Started

### Requirements

- Python 3.x
- pandas
- numpy
- seaborn
- matplotlib
- scipy

### Installation

Clone the repository and install the required dependencies:

```bash
git clone https://github.com/yourusername/gene-expression-analysis.git
cd gene-expression-analysis
pip install -r requirements.txt
