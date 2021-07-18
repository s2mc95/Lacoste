Order 

1)P0 growth
2)P0 settle
3)P14 growth
4)P14 settle 

________________________________________________________________________________________________________________________________________________________________________________________________


CODE DETAILS!!!

1)For P0 growth run 'P0_Growth.m' 
            Line 14-17,37,72,79 tree variables [14-Length ; 15-Width; 16-Height ; 17-tissue height ; 37 - length and width for a smaller area if needed for examination or else same as that of line 14 and 15; 72 - number of root nodes]
            Line 78 learning rate of growth    [eta_wire]  
            Line 80 has stoping variable       [STOP_AT].

            Line 81 : Grow tree wire Function  

            Line 84 : save data line  ; (saved files are named based on STOP_AT: ex 'D_BIG_phase2_just_grown_0.043%_Neurons.mat' has STOP_AT 0.043) 

2)For P0 settle run 'P0Settle.m' 

             Line 2: load the saved data from P0 growth (line 84 of 'P0_Growth.m'  will have the name of the saved mat file).
             Lines 8-12 has the tunable parameters : wire,path,neuron learning rates and path epochs respectively. 
             Line 18: save data line.

3)For P14 growth run 'P14_run.m' 

             Line 3: load the saved data from P0 settle (line 18 of 'P0Settle.m' will have the name of the saved mat file).
             Line 4:  learning rate of growth    [eta_wire]           
           
             Line 5 : Grow tree with neural activity function (grow_tree_with_neural_act_2.m)

            grow_tree_with_neural_act_2 FUNCTION has (4) P14 Settle phase inside it !!!!     
                                 
                       grow_tree_with_neural_act_2 :
                        
                                      Line 5 and 6 : variables to select SOM weights [input_weights & available_weights]  
                                      Line 7 and 8: Height limit to be considered for neural based growth

                                      Line 12 and 232 loads the weight for SOM : For control case 'wt_control_itr'  and for lesioned it is changed to 'wt_lesioned_itr' in both the lines     !!!!
                                      Line 71 Neural Activity Learning Rate

                                      Line 142  Neural Activity Threshold [neural_act_threshold]

                                      Line 230 can be alterd to decide the speed of changing SOM weight parm('wt_control_itr' ) [rem(i,25);i=epoch number ; here every 25th epoch we change]
                                   


                                      Line 250:save data line.


                                      Lines 254,255,256 has the tunable parameters : wire,path and neuron learning rates  respectively 

                                      Lines 258 : path epochs 
                                    
for P14 Settle   [4th PHASE IS INSIDE 3rd PHASE : SETTLE IS INSIDE THE GROWTH!!!!!]
                                
                                      Lines 262 P14 settle function (settle_path_length_witn_neu_act.m) 


                                                              settle_path_length_witn_neu_act :

                                                                               Line 144 has the stopping conditons [(mean(DN2V)<0.15 && mean(DN2V)>0.13) || (RMSE_prev-RMSE_mean(path_epochs,1))<6e-9] .





                       
                                      Lines 267 : Saving Line 


*****in the Final saved file (Line 267 grow_tree_with_neural_act_2 ) Rpos_master will have the values before settle ; Rpos_master_settled will have the values after settle

______________________________________________________________________________________________________________________________________________________________ 


CALCULATION CODES::



* [  In the saved files  Rpos_master is the values before settle ;  Rpos_master_settled are the values after settle]


LSC_VOL_BEFORE_SETTLE.m calculates values before settle   


Line 8 : Load the Control grown Tree data 

Line 63: Load the corresponding Lesioned Tree data 



LSC_VOL_AFTER_SETTLE.m calculates values after settle   


Line 8 : Load the Control grown Tree data 

Line 63: Load the corresponding Lesioned Tree data 





______________________________________________________________________________________________________________________________________________________________
PLOTTING CODE 

portionplot_aftersettle.m 



line 6 : load the saved tree mat file to plot the tree .


