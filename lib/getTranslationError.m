function tranErr = getTranslationError(t_est, t_gt)

    t_est = t_est(:);
    t_gt = t_gt(:);
    tranErr = norm(t_est - t_gt);

end