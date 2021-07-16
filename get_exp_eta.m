function eta = get_exp_eta(epoch, maximum, minimum, max_epoch)
    if epoch < max_epoch
        eta = minimum + (maximum - minimum)*exp(-epoch/(max_epoch));
    else
        eta = minimum;
    end
end