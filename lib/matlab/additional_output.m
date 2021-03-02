function ao = additional_output(additional_output_filename)
ao = containers.Map;

file_contents = fileread(additional_output_filename);
file_lines = strsplit(file_contents, {'\n', '\r'});

for i = 1:length(file_lines)
    thisline = file_lines{i};

    file_data = strsplit(strtrim(thisline));
    mystr = file_data{1};

    if (strcmp(mystr,'!!')),    
        myvar = file_data{2};
		myndims = str2num(file_data{3});
        mydims = zeros(myndims, 1);
        for j=1:myndims
           mydims(j) = str2num(file_data{3+j});
        end
        mydata = zeros(prod(mydims),1);
        for j=1:prod(mydims)
           mydata(j) = str2num(file_data{3+myndims+j}); 
        end
        if myndims > 2,
            mydata = reshape(mydata, mydims); %to be checked...
        end
        ao(myvar) = mydata;
    end
end