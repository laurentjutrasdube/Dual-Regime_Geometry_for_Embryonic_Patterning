# Geometric models for robust encoding of dynamical information into embryonic patterns
Source codes, figure and movies for the paper "Geometric models for robust encoding of dynamical information into embryonic patterns". All codes are written in the `python3` programming language (except for two `Mathematica` notebooks). This repository also contains folders with the source data files, as well as the source codes used to generate the data files.


- **3-gene\_det.ipynb**

    This notebook performs deterministic simulations of the symmetric 3-gene Models 1, 2, 3 and 4. It also performs a bifurcation analysis of these models using the data found in the `XPPAUTO_data` folder, which also contains the `.ode` files used to generate the data with the XPP AUTO software [1]. Figure 2 and Figure 2--figure supplements 1 and 2 show the results obtained with this notebook.
    
    
- **3-gene\_stoch.ipynb**

    This notebook performs stochastic simulations of the symmetric 3-gene Models 1, 2, 3 and 4. It also generates plots of the mutual information using the data found in the `Mutual_info_data` folder, which also contains the `python` codes used to generate the data. Figure 3 and Figure 3--figure supplement 1 show the results obtained with this notebook.
    
    
- **3-gene\_asym.ipynb**

    This notebook performs deterministic simulations of the asymmetric 3-gene Models 1 and 2. It also performs a bifurcation analysis of these models and generates plot of the mutual information using the data found in the `XPPAUTO_data` and `Mutual_info_data` folders, respectively. Figure 5 and Figure 5--figure supplements 1 and 2 show the results obtained with this notebook.
    
    
- **Gene-free\_det.ipynb**

    This notebook performs deterministic simulations of the symmetic gene-free Models 1 and 2. It also performs a bifurcation analysis of these models and generates flow plots using the data found in the `XPPAUTO_data` and `Mathematica_data` folders, respectively. Figure 4, Figure 4--figure supplement 1 and Figure 4--movie supplements 1 and 2 show the results obtained with this notebook.
    
    
- **Gene-free\_stoch.ipynb**

    This notebook performs stochastic simulations of the symmetric gene-free Models 1 and 2. It also generates the mutual information plots using the data found in the `Mutual_info_data` folder. Figure 4--figure supplement 2 shows the results obtained with this notebook.
    
    
    
- **Gene-free\_asym.ipynb**

    This notebook performs deterministic simulations of the asymmetic gene-free Models 1 and 2. It also performs a bifurcation analysis of these models using the data found in the `XPPAUTO_data` folder. Moreover, it generates plots of the flow and of the spatial wave profiles. Figure 6, Figure 6--figure supplement 1, Figure 7 and Figure 7--figure supplement 1 show the results obtained with this notebook.
    
    
- **Hopf\_scenario\_Fig1.ipynb**

    This notebook performs deterministic simulations of the gene network model evolved *in silico* in [2]. Results are shown on Figure 1. It also performs a bifurcation analysis of this model, shown on Figure 1--figure supplement 1.
    
    
- **Infinite-period\_scenario\_Fig1.ipynb**

    This notebook performs deterministic simulations of the infinite-period model of Figure 1, which is a simplified version of the model in the appendix of [3].
    
    
- **Infinite-period\_scenario\_Fig7.ipynb**

    This notebook performs deterministic simulations of the infinite-period model of Figure 7, which is adapted from [4].
    
    
- **Tribolium\_model.ipynb**

    This notebook performs deterministic simulations of the model for Tribolium segmentation from [5]. It also generates flow plots and computes the speed of the cells in phase space. Figure 1--figure supplement 2 shows the results obtained with this notebook.
    

####      

####      

**References**

[1] Bard Ermentrout. Xppaut. In *Computational Systems Neurobiology*, pages 519–531. Springer, 2012.

[2] Paul François, Vincent Hakim, and Eric D Siggia. Deriving structure from evolution: metazoan segmentation. *Molecular systems biology*, 3(1), 2007.

[3] Isabel Palmeirim, Domingos Henrique, David Ish-Horowicz, and Olivier Pourquié. Avian hairy gene expression identifies a molecular clock linked to vertebrate segmentation and somitogenesis. *Cell*, 91(5):639–648, 1997.

[4] Luis G Morelli, Saúl Ares, Leah Herrgen, Christian Schröter, Frank Jülicher, and Andrew C Oates. Delayed coupling theory of vertebrate segmentation. *HFSP journal*, 3(1):55–66, 2009.

[5] Xin Zhu, Heike Rudolf, Lucas Healey, Paul François, Susan J Brown, Martin Klingler, and Ezzat El-Sherif. Speed regulation of genetic cascades allows for evolvability in the body plan specification of insects. *Proceedings of the National Academy of Sciences*, 114(41):E8646–E8655, 2017.

