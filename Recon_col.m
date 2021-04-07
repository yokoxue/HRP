function pm25_recon_sample_c = Recon_col(Sparse_code_c,A,num_b)
rec_Y_l1 = Sparse_code_c;
     tmp_l1 = sort(abs(rec_Y_l1 ),'descend');
     tmp_l1 = tmp_l1(num_b);
     rec_Y_l1 (abs(rec_Y_l1 )<tmp_l1)=0;
     pm25_recon_sample_c = A*rec_Y_l1;
     