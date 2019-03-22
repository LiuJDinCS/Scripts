**简单指南**
  这儿有一系列简单的指南帮助用户开始模拟。更为详细的实例教程可以在这网站找到：[http://www.mdtutorials.com/](http://www.mdtutorials.com/)

# 初学者
  对于那些刚刚开始学习GROMACS或分子动力学模拟的人来说是非常艰难的。强烈建议首先阅读那些为GROMACS提供的各种广泛的文档，以及在感兴趣领域出版的文章。
## 资源
  * GROMACS[参考手册](http://manual.gromacs.org/documentation/2019/reference-manual/index.html)---非常详细的文档，通常也可以作为一个非常好的MD介绍。
  * [流程图](http://manual.gromacs.org/documentation/2019/user-guide/flow.html)---一个蛋白在水何盒子中的典型GROMACS分子动力学流程。
  * 分子动力学模拟和GROMACS介绍([幻灯片](https://extras.csc.fi/chem/courses/gmx2007/Berk_talks/forcef.pdf)，[视频](http://tv.funet.fi/medar/showRecordingInfo.do?id=/metadata/fi/csc/courses/gromacs_workshop_2007/IntroductiontoMolecularSimulationandGromacs_1.xml))---力场，积分以及温度和压力的控制（Berk Hess）

# 添加残基到力场中

## 添加新的残基
  如果你有需要将新的残基添加到已经存在的力场中去以便能够使用[pdb2gmx](http://manual.gromacs.org/documentation/2019/onlinehelp/gmx-pdb2gmx.html#gmx-pdb2gmx)，或者修改存在的残基，那么有几个文件你应该修改。你必须参阅[参考手册](http://manual.gromacs.org/documentation/2019/reference-manual/index.html)对于所需格式的描述部分。执行以下步骤：
  * 向你选择的力场中的[rtp](http://manual.gromacs.org/documentation/2019/reference-manual/file-formats.html#rtp)文件中添加残基.你也可以复制一个已经存在的残基，重命名和适当地修改它，或者你可能需要使用额外的拓扑生成工具并调整其为[rtp](http://manual.gromacs.org/documentation/2019/reference-manual/file-formats.html#rtp)格式。

## 修改一个力场

# 水溶剂

# 非水溶剂

## 制作非水溶剂盒子

# 混合溶剂

# 制作二硫键

# GROMACS运行膜的模拟

## 运行膜模拟

## 用genbox添加水

## 扩充材料

# 新分子参数化

## 外来物种

# 平均力势

# 单点能

# 碳纳米管

## Robert Johnson的技巧
  采纳了Robert Johnson在gmx-users mailing list上的帖子。
  * 要绝对保证拓扑文件中“末端”碳原子是共享一个键。
  * [gmx grompp](http://manual.gromacs.org/documentation/2019/onlinehelp/gmx-grompp.html#gmx-grompp)时的[mdp](http://manual.gromacs.org/documentation/2019/reference-manual/file-formats.html#mdp)输入文件要使用**periodic_molecules = yes**。
  * 即使拓扑文件正确，如果你把碳纳米管放置在错误尺寸的盒子里将使得扭曲发生。因此使用[VMD](http://www.ks.uiuc.edu/Research/vmd/)可视化纳米管和它的周期性成像，确保镜像之间的空隙正确。如果空隙太小或太大，在纳米管上将有巨大的压力导致其扭曲或拉伸。
  * 沿着碳纳米管的轴线方向上不要使用压力耦合。事实上，为了调试目的，如果发生某些错误，可能先最好关掉整个压力耦合直到你弄清楚是什么问题。
  * 当使用具有特定力场的[x2top](http://manual.gromacs.org/documentation/2019/onlinehelp/gmx-x2top.html#gmx-x2top)时，假设分子的连接性。如果是周期性的，亦或者是非周期性上面有氢原子，那么你的纳米管末端的碳原子最多只能和两个其他分子相结合。
  * 你可以使用[x2top](http://manual.gromacs.org/documentation/2019/onlinehelp/gmx-x2top.html#gmx-x2top)工具中的 **-pbc** 选项产生一个无限长的纳米管，这里[x2top](http://manual.gromacs.org/documentation/2019/onlinehelp/gmx-x2top.html#gmx-x2top)将识别为末端C原子共享一个化学键。因此，当你使用[grompp](http://manual.gromacs.org/documentation/2019/onlinehelp/gmx-grompp.html#gmx-grompp)时不会得到关于单个碳原子的报错。

## Andrea Minoia的指南

  GROMACS模拟碳纳米管（也可以参看：[http://www.webcitation.org/66u2xJJ3O](http://www.webcitation.org/66u2xJJ3O)）

  **原文翻译**
  **使用GROMACS模拟碳纳米管**
  网上有许多关于GROMACS和CNTs的文献（[1](http://www.gromacs.org/Documentation/How-tos/Carbon_Nanotube),[2](http://cs86.com/CNSE/SWNT.htm),[3](http://machine-phase.blogspot.com/2009/04/single-wall-carbon-nanotubes-in-403.html)，[...](http://www.google.be/search?hl=en&source=hp&q=cnt+gromacs&btnG=Google+Search&meta=&aq=f&oq=)）。即使读了很多，我也没有找到我想要的：一个清楚、简单的模拟单壁碳纳米管和多壁碳纳米管（SWNT，MWNT）指南，周期或者非周期性。
  我花了一段时间才找到步骤做这些我想做的事情：构建一个周期性CNT并创建GROMACS的拓扑文件。

  ## 构建和准备CNT
  第一步很明显是构建CNT，列入采用 [YASC buildCstruct](http://chembytes.wikidot.com/buildcstruct)。使用这个插件可以构建SWCNT和MSWCNT，包括armchair或zigzag型，有限的或周期性的。从 v1.1版本开始，**buildCstruct**能够直接保存为grmacs的GRO格式文件。
  **注意要有一个足够大的晶胞容纳所有的CNT，否则你将不能产生拓扑文件。你能够使用ngmx，molden或者vnd的pbc检查结构。而且，也要注意选择一个足够大的盒子去容纳其他你想和CNT接触的分子，例如聚合物链**
  
  ## 为x2top准备的文件
  x2top需要一堆文件来创建拓扑文件。我想使用oplsaa力场。但是我不想搞乱**share/gromacs/top**目录中的原始ffoplsaa文件，因此我自己创建了我需要的文件

  ### .n2t itp
  此文件用于在拓扑中转换名称。根据它的连接性读取结构中原子的名称（例如成键数，键长和原子结合），根据.n2t文件中定义的尝试去猜合适的原子类型。我的基于oplsaa力场的石墨烯和CNT的.n2t文件如下：
```
; Oplsaa-based n2t for carbon-based structures such as CNTs and graphenes
; Andrea Minoia
H    HJ    0.00       1.008  1    C 0.109                              ;Hydrogen
C    CJ    0.00      12.011  3    C 0.142   H 0.109   H 0.109 ;Periferic C
C    CJ    0.00      12.011  3    C 0.142   C 0.142   H 0.108 ;Periferic C
C    CJ    0.00      12.011  1    C 0.142                          ;Internal/periodic C
C    CJ    0.00      12.011  2    C 0.142   C 0.142            ;Internal/periodic C
C    CJ    0.00      12.011  3    C 0.142   C 0.142   C 0.142 ;Internal/periodic C
```
  我决定把所用的电荷设置为0，但是你可能想不一样。

  ### .itp文件
  我更喜欢不让x2top将力场参数放入拓扑文件中，我已经定义了我自己的.itp文件：
```
; Oplsaa-based force field for carbon-based structures such as CNTs and graphenes
; Andrea Minoia

[ defaults ]
; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ
1               3               yes             0.5     0.5
; parameters are taken from the OPLS force field

[ atomtypes ]
; The charges here will be overwritten by those in the rtp file
; name       mass      charge    ptype      sigma      eps
  CJ   6     12.01100     0.000       A    3.55000e-01  2.92880e-01 ;opls_147 naftalene fusion C9
  HJ   1      1.00800     0.000       A    2.42000e-01  1.25520e-01 ;opls_146 HA hydrogen benzene. I have set the charges zero

[ bondtypes ]
; i    j func        b0          kb
  CJ    CJ      1    0.14000   392459.2   ; TRP,TYR,PHE
  CJ    HJ      1    0.10800   307105.6   ; PHE, etc.

[ angletypes ]
  CJ     CJ     CJ      1   120.000    527.184   ; PHE(OL)
  CJ     CJ     HJ      1   120.000    292.880   ;
  HJ     CJ     HJ      1   117.000    292.880   ; wlj from HC-CM-HC

[ dihedraltypes ]
  CJ     CJ     CJ     CJ      3     30.33400   0.00000 -30.33400   0.00000   0.00000   0.00000 ; aromatic ring
  HJ     CJ     CJ     HJ      3     30.33400   0.00000 -30.33400   0.00000   0.00000   0.00000 ; aromatic ring
  HJ     CJ     CJ     CJ      3     30.33400   0.00000 -30.33400   0.00000   0.00000   0.00000 ; aromatic ring
```
目前没有二面角。
  ### .rtp 文件
  最后，我需要创建oplsaa的残基拓扑文件，这里还没有定义任何残基。
```
; New format introduced in Gromacs 3.1.4.
; Dont use this forcefield with earlier versions.

; Oplsaa-based rtp for carbon-based structures such as CNTs and graphenes
; Andrea Minoia

; NB: OPLS chargegroups are not strictly neutral, since we mainly
; use them to optimize the neighborsearching. For accurate simulations
; you should use PME.

[ bondedtypes ]
; Col 1: Type of bond
; Col 2: Type of angles
; Col 3: Type of proper dihedrals
; Col 4: Type of improper dihedrals
; Col 5: Generate all dihedrals if 1, only heavy atoms of 0.
; Col 6: Number of excluded neighbors for nonbonded interactions
; Col 7: Generate 1,4 interactions between pairs of hydrogens if 1
; Col 8: Remove propers over the same bond as an improper if it is 1
; bonds  angles  dihedrals  impropers all_dihedrals nrexcl HH14 RemoveDih
     1       1          3          1        0         3      1     0
```
  ### FF.dat文件
  这里我指定了力场的名字：
```
1
ffcnt_oplsaa
```
  ## 使用x2top创建拓扑
  一旦所有文件准备好，就可以使用x2top产生拓扑：
```
x2top -f file.gro -o topol.top -ff cnt_oplsaa -name CNT -noparam -pbc
```
  **-pbc**选项将产生那些周期性边界条件的所有键，角度和二面角（这就是为什么.gro文件中需要指定合适大小的盒子单元）。如果你对纳米管不是周期性的，那就不需要 **-pbc** 。
  ## 矫正拓扑
  最后一步是将用来描述二面体的函数从**1**更改为**3**。
  ## .mdp文件
  在你的mdp文件中打开周期性选项：
```
periodic_molecules=yes
pbc = xyz
```
  这就是全部...现在你已经可以模拟周期性碳纳米管或碳基结构。
  ## 溶解CNT
  为了溶解CNT，我们能够使用GROMACS的genbox工具（5.0版本以上改为了solvate）采用=cp（纳米管的gro文件）和-cs（溶解的盒子gro文件）。
  通常，我首先在NPT下平衡溶剂盒子，然后创建了一个比CNT盒子更大的超胞。这是因为溶质盒子（-cp）将成为最终系统的盒子。
  genbox将尝试用溶剂分子填充盒子空隙，甚至在CNT的内部。如果这不可取，我们能写一个创建索引文件的脚本，该文件不包含CNT内部分子中原子的原子索引。然后，使用editconf去产生最终体系：
```
editconf -f starting_system.gro -n index.ndx -o final_system.gro
```
  ## 使用YSAC nosolincnt去除CNT内部的原子
  与editconf一起使用的引索文件是为了去除CNT内部的溶剂分子，这能够佷容易通过YSAC nosolincnt编写，它是一个VMD插件。
  用VMD读取体系后，打开TCL控制台，执行插件：
```
source /path/of/the/script/nosolincnt.tcl
```
  将提示你插入CNT的片段数和CNT半径（埃米）。这些信息将用于创建一个出圆柱形，根据CNT质心X，Z坐标得到的半径。当视图窗口处于激活状态时，你能按 **C** 去检查选择的结果。如果一切正确，然后按 **W** 将会写入引索文件而不需要选择原子。
  在下面的图片中，分别显示了溶剂化后的碳纳米管、要除去的分子和最终的体系
  {% asset_img cnt.jpg CNT体系 %}

-----
  包括了使用OPLS-AA力场设置一个CNT模拟的详细步骤。一个CNT的结构可以很容易通过[buildCstruct](http://chembytes.wikidot.com/buildcstruct)（此Python脚本也可以加氢）或者通过[TubeGen Online](http://turin.nss.udel.edu/research/tubegenonline.html)（只需要复制和粘贴PDB输出到文件并命名为cnt.pdb）。
  为了能够用GROMACS模拟，你可能要执行以下操作：
  * 新建一个cnt_oplsaa.ff目录文件夹
  * 在此目录下，创建如下文件，以下步骤来源于上面的指南：
    * forcefield.itp来自[itp](http://manual.gromacs.org/documentation/2019/reference-manual/file-formats.html#itp)部分文件。
    * atomnames2types.n2t来自[n2t](http://manual.gromacs.org/documentation/2019/reference-manual/file-formats.html#n2t)部分文件。
    * aminoacids.rtp来自[rtp](http://manual.gromacs.org/documentation/2019/reference-manual/file-formats.html#rtp)部分文件。
  * 使用传统的力场产生拓扑（cnt_oplsaa.ff目录文件必须与[gmx x2top](http://manual.gromacs.org/documentation/2019/onlinehelp/gmx-x2top.html#gmx-x2top)命令执行的目录相同，亦或者它能够在GMXLIB路径中找到），采用[gmx x2top](http://manual.gromacs.org/documentation/2019/onlinehelp/gmx-x2top.html#gmx-x2top)的 **-noparam** 选项，不使用来自命令行（-kb, -ka, -kd）指定的bond/angle/dihedral等常数，而是依赖力场文件；然而，这就需要下一步（确定二面角函数）
```
gmx x2top -f cnt.gro -o cnt.top -ff cnt_oplsaa -name CNT -noparam
```
二面角的函数类型通过[gmx x2top](http://manual.gromacs.org/documentation/2019/onlinehelp/gmx-x2top.html#gmx-x2top)设置为“1”，但是力场文件指定的函数类型为“3”。因此，替换拓扑文件中 **[ dihedrals ]** 部分的”1“为”3“。一种快速的方式是使用 **sed**（但是你核能需要调整你的操作系统；也要手动看一下top文件核查一下你更改的二面角函数类型）：
```
sed -i~ '/\[ dihedrals \]/,/\[ system \]/s/1 *$/3/' cnt.top
```
一旦你有了拓扑文件，你就可以设置你的系统了。例如，一个简单的 真空模拟（在[em.mdp](http://manual.gromacs.org/documentation/2019/reference-manual/file-formats.html#mdp)和[md.mdp](http://manual.gromacs.org/documentation/2019/reference-manual/file-formats.html#mdp)中使用你自己的参数）：
放入一个稍微大的盒子中：
```
gmx editconf -f cnt.gro -o boxed.gro -bt dodecahedron -d 1
```
真空中的能量最小化：
```
gmx grompp -f em.mdp -c boxed.gro -p cnt.top -o em.tpr
gmx mdrun -v -deffnm em
```
真空中的MD：
```
gmx grompp -f md.mdp -c em.gro -p cnt.top -o md.tpr
gmx mdrun -v -deffnm md
```
查看轨迹：
```
gmx trjconv -f md.xtc -s md.tpr -o md_centered.xtc -pbc mol -center
gmx trjconv -s md.tpr -f md_centered.xtc -o md_fit.xtc -fit rot+trans
vmd em.gro md_fit.xtc
```

# 可视化软件

## 拓扑成键 vs 渲染的键

# 提取轨迹信息

# 外部工具执行轨迹分析

# 绘制数据

## 软件

# 胶束聚类

## 修改力场