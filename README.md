# InverseSolution_for_ERP

## Inverse solution 
Structural processing and source time series are computed using CiftiStorm pipeline (Areces-Gonzalez et al., 2024), a MATLAB toolbox for spectral analysis and inverse solution of EEG/MEG data. 
Spectral analysis is via K-order orthogonal wavelet decomposition of the data (Eq.1) employing maximal overlap discrete wavelet transform method (MODWT) (Percival et al. 1997).

〖MODWT〗_K: v(t)→v(t,k);t∈[0,..,T],k∈[1,..,K+1] 				(Eq.1)

In our convention v(t) is a vector function in the discrete time domain t and v(t,k) (wavelet transform) a vector function in the discrete time domain t and wavelet domain k. Wavelets are waveforms peaking around bands of frequencies described by the rule (Fs/2^(k+1) ,├ Fs/2^k ]┤ that involves the wavelet order k and the sampling rate Fs. 

In our settings, a sampling rate Fs=128Hz, and with a maximum wavelet order of K=7, waveforms peak around the fundamental EEG frequency bands: (k=1 ) Gamma 32-64Hz, (k=2) Beta 16-32Hz, (k=3) Alpha 8-16Hz, (k=4) Theta 4-8Hz, (k=5) High Delta 2-4Hz, (k=6) Mid Delta 1-2Hz, (k=7)  Low Delta 0.5-1Hz, and the low frequency drift  (k=8) <0.5Hz.

## Inverse solutions are then obtained via the Spectral Structured Sparse Bayesian Learning method (SSSBL) (Paz-Linares et al., 2023). 
The SSSBL inverse solution (Eq.2) is via an optimal quasilinear operator applied to the wavelet transform data v(t,k) that obtains ι(t,k). 
Note that this approach optimizes an operator to target the source activity patterns produced by cortical oscillatory networks at each of the major bands of frequencies in resting state EEG.   

SSSBL H_ιv (S_vv (k)): v(t,k)→ ι(t,k);t∈[0,..,T] 					(Eq.2)

where S_vv (k)=(∑_t▒〖v^† (t,k)v(t,k) 〗)⁄T covariance of the sensor waveforms in time.

## Here, the quasilinear operator H_ιv (S_vv (k)) denotes a nonlinear matrix function of a matrix argument, the wavelet data covariance S_vv (k). 

Quasilinear is a quality of inverse solutions like SSSBL, or others that pursue strong regularization via nonlinear priors, controlling spatial distortions or leakage of the source time series (Gonzalez-Moreira, 2018), while avoiding nonlinear distortions in time of frequency domains. 

### Calculations of the SSSBL operator (Eq.3) optimize the following Bayesian cost function.

H_ιv (S_vv (k))=〖argmax〗_(H_ιv ) {Q(S_vv (k),H_ιv )}						(Eq.3)

where Q(S_vv (k),H_ιv )=R(S_vv (k),H_ιv )+P(S_vv (k),H_ιv ) cost function (log-posterior), with:

R(S_vv (k),H_ιv )= tr((I_vv-L_vι H_ιv ) S_vv (k) (I_vv-L_vι H_ιv )^† ) log-likelihood

 P(S_vv (k),H_ιv )= λ_1 〖tr〗^(1/2) (H_ιv S_vv (k) H_ιv^† )+λ_2 tr(H_ιv S_vv (k) H_ιv^† ) log-prior

### The choice of a strong nonlinear prior model above is the denominated Elastic Net Nuclear Norm, which targets S_ιι (k)  the covariance of source waveforms in time.

where S_ιι (k)=(∑_t▒〖ι^† (t,k)ι(t,k) 〗)⁄T=H_ιv S_vv (k) H_ιv^†

### Finally, the inverse solution (source time series ι(t)) summarizes in (Eq.4), applying the inverse MODWT to the source waveforms ι(t,k).

〖IMODWT〗_K: ι(t,k)→ι(t);t∈[0,..,T],k∈[1,..,K+1] 					(Eq.4)

## References 

1- Areces-Gonzalez, A., Paz-Linares, D., Riaz, U., Wang, Y., Li, M., Razzaq, F.A., Bosch-Bayard, J.F., Gonzalez-Moreira, E., Lifespan Brain Chart Consortium (LBCC), Global Brain Consortium (GBC) and Cuban Human Brain Mapping Project (CHBMP), 2024. Ciftistorm pipeline: facilitating reproducible eeg/meg source connectomics. Frontiers in Neuroscience, 18, p.1237245.
2- Percival, D.B. and Walden, A.T., 2000. Wavelet methods for time series analysis (Vol. 4). Cambridge university press.(Percival et al. 1997)
3- Paz-Linares, D., Gonzalez-Moreira, E., Areces-Gonzalez, A., Wang, Y., Li, M., Vega-Hernandez, M., Wang, Q., Bosch-Bayard, J., Bringas-Vega, M.L., Martinez-Montes, E. and Valdes-Sosa, M.J., 2023. Minimizing the distortions in electrophysiological source imaging of cortical oscillatory activity via Spectral Structured Sparse Bayesian Learning. Frontiers in Neuroscience, 17, p.978527.
4- Gonzalez-Moreira, E., Paz-Linares, D., Areces-Gonzalez, A., Wang, R., Bosch-Bayard, J., Bringas-Vega, M.L. and Valdes-Sosa, P.A., 2018. Caulking the leakage effect in MEEG source connectivity analysis. arXiv preprint arXiv:1810.00786.

## Dependencies
1- EEGLAB (version 2024 or highier)

## Test data
1- Download the test data from [[Here]](https://lstneuro-my.sharepoint.com/:f:/g/personal/ariosky_neuroinformatics-collaboratory_org/IgAW2KHcsUCKT5L1pXGg5tSFAbRyuXSyNBcjEuSt2EUATTw?e=k6sO5x)
2- Unzip the file in the project directory
 

## Params configuration
### General params file
- app/general_params.json  (for the text data keep the general params configuration )
- Define the path to EEGLAB
"dependencies":{
        "eeglab":{
            "base_path":"/root/path/to/EEGLAB"
        }
    }

### Event params file
app/event_params.json
- select tag (make true the events or conditions you wantto analize ex ("get":true))
- slection tag (select an option)
    - "type": "segments" ( it will divide the data in segment of the same condition)
    - "type": "trials" (it  will select the trials based on nTrials and time 
    - "type": "continue" ( it will select the segments of the same condition and join them as continues data)