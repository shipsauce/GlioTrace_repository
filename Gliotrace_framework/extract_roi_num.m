function roi_num = extract_roi_num(str)
    match = regexp(str, 'roi_(\d+)', 'tokens');
    roi_num = str2double(match{1}{1});
end