# HFOApp

A MATLAB-based graphic user interface for high frequency oscillations marking.

Developed using APP Designer (MATLAB R2020b).

Related to Guangyu Zhou, Torben Noto, Arjun Sharma, Qiaohan Yang, Karina A. González Otárula, Matthew Tate, Jessica W. Templer, Gregory Lane and Christina Zelano. HFOApp: A MATLAB graphical user interface for high frequency oscillation marking. eNeuro 20 September 2021, ENEURO.0509-20.2021; DOI: https://doi.org/10.1523/ENEURO.0509-20.2021 


Comments and questions: guangyu.zhou@northwestern.edu




05/01/2023. Fix bug in HFOAutoDetectorHil.m       
          This bug affects the Hilbert Detector that was called by 
          Options -> Automatic HFO detection -> Hilbert Detector (HFOApp), only when "Epoch (s)" was edited.


03/06/2024. An error occurred when select to display events on current channel. 
	The bug was fixed by changing line 2761: if app.ShowEventOnAllChannelButton.Value == 0 && ~ismember( cur_label, app.uidata.current_channel) 
	to         
	if app.ShowEventOnAllChannelButton.Value == 0 && ...
                    (isempty( app.uidata.current_channel) || ~ismember( cur_label, app.uidata.current_channel( :, 1)))

      HFOApp.mlapp was updated.

