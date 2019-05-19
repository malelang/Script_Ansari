function DL_model_builder(alarm)

	addpath(genpath('./deep_learning_toolbox/'));

    no_folds_test = 5;
    no_folds_val = 5;
    
    train_func = @DL_train;
    test_func = @DL_test;
    
    params = {50:50:300};

    load([alarm '_wrapped_data_new_norm.mat']);
    
    signals = wrapped_data.signal;
    signals(signals(:,end)==1,end) = 0;
    signals(signals(:,end)==2,end) = 1;    
    
    signals = [signals ; repmat(signals(signals(:,end)==1,:),6,1)];
    
    data = signals(:,[1:2:250, end]);
    
    data0 = data(data(:,end)==0,:);
    data1 = data(data(:,end)==1,:);
    
    rperm0 = randperm(length(data0));
    rperm1 = randperm(length(data1));
    data0 = data0(rperm0,:);
    data1 = data1(rperm1,:);
    
    fold_size0 = ceil(size(data0,1)/no_folds_val);
    fold_size1 = ceil(size(data1,1)/no_folds_val);    
    
    accs = [];
    accs0 = [];
    accs1 = [];
    tk = tic;
    parfor fi = 1:no_folds_val
        fold_test0 = data0(fold_size0*(fi-1)+1:min(fold_size0*fi,size(data0,1)),:);
        fold_test1 = data1(fold_size1*(fi-1)+1:min(fold_size1*fi,size(data1,1)),:);

        fold_train0 = data0([1:fold_size0*(fi-1),min(fold_size0*fi,size(data0,1))+1:end],:);
        fold_train1 = data1([1:fold_size1*(fi-1),min(fold_size1*fi,size(data1,1))+1:end],:);            

        fold_test = [fold_test0 ; fold_test1];
        fold_test = fold_test(randperm(size(fold_test,1)),:);

        fold_train = [fold_train0 ; fold_train1];
        fold_train = fold_train(randperm(size(fold_train,1)),:);            

        best_model = crossvalidation(train_func, test_func, fold_train, params, no_folds_test);

        [acc, acc0, acc1] = test_func(best_model,fold_test);
        
        accs = [accs acc];
        accs0 = [accs0 acc0];
        accs1 = [accs1 acc1];
    end    
    toc(tk);
    acc = sum(accs)/no_folds_val;
    acc0 = sum(accs0)/no_folds_val;
    acc1 = sum(accs1)/no_folds_val;
    
    [best_model, best_pars] = crossvalidation(train_func, test_func, data(randperm(size(data,1)),:), params, no_folds_test);    
    
    disp(best_pars);
    disp(acc);
    disp(acc0);
    disp(acc1);
    
    save([alarm '_DL_model'],'best_model');
end


function [best_model, best_pars] = crossvalidation(train_func, test_func, data,params,no_folds)

    data0 = data(data(:,end)==0,:);
    data1 = data(data(:,end)==1,:);
    
    num_total = num_pars(params);
    pars = zeros(num_total,length(params));
    for pi = 1:length(params)
        num = num_pars(params(pi:end));
        num2 = num_pars(params(pi+1:end));
        
        for pj = 1:num_total/num
            for pk = 1:length(params{pi})
                pars((pj-1)*num+(pk-1)*num2+1:(pj-1)*num+pk*num2,pi) = params{pi}(pk);
            end
        end
    end
    
    rperm0 = randperm(length(data0));
    rperm1 = randperm(length(data1));
    data0 = data0(rperm0,:);
    data1 = data1(rperm1,:);
    
    fold_size0 = ceil(size(data0,1)/no_folds);
    fold_size1 = ceil(size(data1,1)/no_folds);

    results = zeros(size(pars,1),size(pars,2)+1);
    parfor pi = 1:size(pars,1)
        acc = 0;
        for fi = 1:no_folds
            fold_test0 = data0(fold_size0*(fi-1)+1:min(fold_size0*fi,size(data0,1)),:);
            fold_test1 = data1(fold_size1*(fi-1)+1:min(fold_size1*fi,size(data1,1)),:);

            fold_train0 = data0([1:fold_size0*(fi-1),min(fold_size0*fi,size(data0,1))+1:end],:);
            fold_train1 = data1([1:fold_size1*(fi-1),min(fold_size1*fi,size(data1,1))+1:end],:);            
            
            fold_test = [fold_test0 ; fold_test1];
            fold_test = fold_test(randperm(size(fold_test,1)),:);
            
            fold_train = [fold_train0 ; fold_train1];
            fold_train = fold_train(randperm(size(fold_train,1)),:);            
            
            model = train_func(fold_train,pars(pi,:));
            
            acc = acc + test_func(model,fold_test);            
        end
        acc = acc / no_folds;
        results(pi,:) = [pars(pi,:), acc];
    end

    [~,idx] = max(results(:,end));
    best_pars = results(idx,1:end-1);
    
    data = data(randperm(size(data,1)),:);
    best_model = train_func(data,best_pars);
end

function par_no = num_pars(pars)
    par_no = 1;
    for i = 1:length(pars)
        par_no = par_no*length(pars{i});
    end
end


function nn = DL_train(data, pars)
    opts.numepochs =  100;
    opts.batchsize = 100;

    data = data(1:floor(size(data,1)/opts.batchsize)*opts.batchsize,:);
    
    nn = nnsetup([size(data,2)-1, pars, 2]);

%    tk = tic;
    [nn, L] = nntrain(nn, data(:,1:end-1), [data(:,end), 1-data(:,end)], opts);    
%    toc(tk);
end

function [acc, acc0, acc1] = DL_test(nn,data)
    labels = nnpredict(nn, data(:,1:end-1));
    labels(labels==2) = 0;
    
    acc = 1 - sum(abs(labels - data(:,end)))/length(labels);
    acc0 = 1 - sum(abs(labels(data(:,end)==0) - 0))/sum(data(:,end)==0);
    acc1 = 1 - sum(abs(labels(data(:,end)==1) - 1))/sum(data(:,end)==1);
end