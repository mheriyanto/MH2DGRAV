![GitHub language count](https://img.shields.io/github/languages/count/mheriyanto/MH2DGRAV)
![GitHub top language](https://img.shields.io/github/languages/top/mheriyanto/MH2DGRAV)
![GitHub repo size](https://img.shields.io/github/repo-size/mheriyanto/MH2DGRAV)
![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/mheriyanto/MH2DGRAV)
![GitHub last commit](https://img.shields.io/github/last-commit/mheriyanto/MH2DGRAV.svg)
[![HitCount](http://hits.dwyl.com/mheriyanto/MH1DDC.svg)](http://hits.dwyl.com/mheriyanto/MH2DGRAV)
[![LinkedIn](https://img.shields.io/badge/-LinkedIn-black.svg?style=flat&logo=linkedin&colorB=555)](https://id.linkedin.com/in/mheriyanto)

# MH2DGRAV
MH2DGRAV is continuous two-dimension inversion of Gravity data based on Talwani formulation using very fast simulated annealing (VFSA) in MATLAB.

These were scripts that were used to implement our paper: **W. Srigutomo, M. Heriyanto, and M. Hilmi Aufa. Gravity Inversion of Talwani Model using Very Fast Simulated Annealing. Journal of Mathematical and Fundamental Sciences, Vol. 51, No. 2, 2019, 177-190. doi: 10.5614/j.math.fund.sci.2019.51.2.7** ([**PDF**](http://journals.itb.ac.id/index.php/jmfs/article/download/8270/4113)). I presented this paper on International Conference on Mathematics and Natural Sciences 2016. Nov 2, 2016 ([**SLIDE**](https://figshare.com/articles/Gravity_Inversion_Based_on_Talwani_Formulation_Using_Very_Fast_Simulated_Annealing_Schemes/5296225)). These scripts contain two main scripts: [forward](https://github.com/mheriyanto/MH2DGRAV/tree/master/forward) and [VFSA inversion](https://github.com/mheriyanto/MH2DGRAV/tree/master/vfsa_inversion): horizontal model. Tutorial video on **YouTube**: [VFSA](https://youtu.be/4UbVILIUDwI). 

I hope these scripts can help students to enter research on geophysical inversion. Any updates about these scripts can be seen in [my blog](https://mheriyanto.wordpress.com/mh2dgrav/): https://mheriyanto.wordpress.com/mh2dgrav.

<ins>**Forward result**</ins>

+ **Kernel Matrix for forward modeling**

<p align="center">
<img src="https://github.com/mheriyanto/MH2DGRAV/blob/master/forward/noise%2010%25/horizontal%20model/Matrix%20A%20Plot.png" width="70%">
</p>

+ **Data with 10% noise by horizontal model**

<p align="center">
<img src="https://github.com/mheriyanto/MH2DGRAV/blob/master/forward/noise%2010%25/horizontal%20model/Plot%20forward%20of%20horizontal%20model.png" width="70%">
</p>

+ **Data with 10% noise by vertical model**

<p align="center">
<img src="https://github.com/mheriyanto/MH2DGRAV/blob/master/forward/noise%2010%25/vertical%20model/Plot%20forward%20of%20vertical%20model.png" width="70%">
</p>

<ins>**VFSA Inversion result**</ins>

+ **Data with 10% noise by horizontal model**

<p align="center">
<img src="https://github.com/mheriyanto/MH2DGRAV/blob/master/results/horizontal_model/Plotting%20Final%20VFSA.png" width="50%">
<img src="https://github.com/mheriyanto/MH2DGRAV/blob/master/results/horizontal_model/RMSE2Iteration.png" width="30%">
</p>

+ **Data with 10% noise by vertical model**

<p align="center">
<img src="https://github.com/mheriyanto/MH2DGRAV/blob/master/results/vertical_model/Plotting%20Final%20VFSA.png" width="50%">
<img src="https://github.com/mheriyanto/MH2DGRAV/blob/master/results/vertical_model/RMSE2Iteration.png" width="30%">
</p>

## Usage

```console
$ git clone https://github.com/mheriyanto/MH2DGRAV.git
$ cd MH2GRAV
$ cd vfsa_inversion
$ octave Grav2DInvVFSAHorisontal.m
```

## License
MH2DGRAV is released under the MIT License (refer to the [LICENSE](https://github.com/mheriyanto/MH2DGRAV/blob/master/LICENSE) file for details).

## Citation
If you find this project useful for your research, please use the following BibTeX entry.

```BibTeX
@article{srigutomo2019gravity,
    title={Gravity Inversion of Talwani Model using Very Fast Simulated Annealing},
    author={Srigutomo, Wahyu and Heriyanto, Mohammad and Aufa, Muhamad Hilmi},
    journal={Journal of Mathematical and Fundamental Sciences},
    volume={51},
    number={2},
    pages={177--190},
    year={2019}
}
```
