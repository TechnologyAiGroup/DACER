# The DACER Project

This is the website for the code and example data in the paper, "DACER: Diagonal-Adaptive CMA-ES with Early-Exit for Reliability-based Attacks on Strong PUFs".

Now it contains the attack method for multiple PUF structures in the folders named after the type of PUF.

## Usage

To start the attack, please go to the corresponding PUF folder and excute

```
matlab main.m
```

The attack program will automatically generate CRPs and start attacking, if you want to use the example dataset, uncomment the lines about reading datasets from the file.

In the `main.m` folder, you can freely define the following parameters to change the PUF structure and attack method.
* The bit number `n`, structure size `nXOR` and noise level `sigmaNoise` of the PUF.

* The training data size `nTrS` and the test data size `nTeS` used for the attack.

* `flag_diag` indicates whether only the elements on the diagonal are used for updates in the CMA-ES algorithm. Enabling this parameter can reduce the running time of each round of the CMA-ES algorithm. Setting it to 0 indicates off and 1 indicates on.

* `flag_early` indicates whether the CMA-ES algorithm is prematurely ended when the variation range of the fitness values of the optimal candidate solutions in ten consecutive rounds is less than 10^-6. Enabling this parameter when the training set is sufficient can reduce the running time of each round of the CMA-ES algorithm. Setting it to 0 indicates off and 1 indicates on.

* `flag_fitloss` indicates whether to add a penalty term to limit the value of the optimal candidate solution in each round of the CMA-ES algorithm to be less than the given value `dis_mu`.  This parameter needs to be used in combination with `loss_func`, `loss_alpha`, `dis_mu` and `dis_sig`.  Under normal circumstances, `loss_alpha` represents the coefficient of the penalty term, and `dis_mu` and `dis_sig` are related to the specific PUF structure.  Enabling `flag_fitloss` can prevent the CMA-ES algorithm from converging to the wrong solution and reduce the number of running rounds of the algorithm.  Setting it to 0 indicates off and 1 indicates on

In specific attacks, `flag_diag`, `flag_early` and `flag_fitloss` can be freely combined and used, which will produce different attack effects. For the best attack effect, we recommend enabling all of them.

The `{PUF_type}_ATTACK.m` file under each folder provides the specific attack steps for different PUF structures. The default settings can be used directly without any changes.

The `cmaes.m` file provides the specific code of the optimized and improved CMA-ES algorithm. The default settings can be used directly without any changes.

Other common functions used during the attack process are placed in the `./puf_util` folder.

For Transfer_DACER, please go to the corresponding folder and excute

```
matlab Transfer_DACER.m
```
The main modifiable parameters include The bit number 'n', structure size 'nXOR', noise level 'sigmaNoise' of the PUF, the source domain training data size `nTrS` and the target domain training data size `nTrS_t`.

Due to the randomness of the CMA-ES algorithm, please run the above attack code multiple times to ensure success.

For more details, please check the specific MATLAB file.

