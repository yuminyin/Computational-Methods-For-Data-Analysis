clear all; close all; clc

[images, labels] = mnist_parse('train-images-idx3-ubyte', 'train-labels-idx1-ubyte');
[test_images, test_labels] = mnist_parse('t10k-images-idx3-ubyte', 't10k-labels-idx1-ubyte');

%% reshape each image and perform an SVD analysis
[~,~,num_of_images] = size(images);
images_reshape = reshapeimages(images);

% Subtract the mean of each row from the data 
mean_of_row = mean(images_reshape, 2); 
images_demean = images_reshape - repmat(mean_of_row, 1, num_of_images);

% demean the test data
test_images_reshape = reshapeimages(test_images);
[~,~,num_of_test_images] = size(test_images);
test_images_reshape_demean = test_images_reshape - repmat(mean_of_row, 1, num_of_test_images); % demean the test data

% perform an SVD analysis
[U,S,V] = svd(images_demean,'econ');

%% Plot the singular value spectrum

plot(diag(S),'ko','Linewidth',2)
title('Singular Value Spectrum')
set(gca,'Fontsize',14,'Xlim',[0 784])
xlabel('singular indices')
ylabel('singular values')

%% Project onto three selected V-modes 

Projection = U(:,[1,2,3])' * images_reshape;
for i = 0:9
    Projection_digits = Projection(:, find(labels == i));
    plot3(Projection_digits(1,:), Projection_digits(2,:), Projection_digits(3,:), 'o'); hold on
end
legend('0', '1', '2', '3', '4', '5', '6', '7', '8', '9')
title('Projection onto Three Selected V-modes')
set(gca, 'fontsize', 14)
xlabel('1st mode')
ylabel('2nd mode')
zlabel('3rd mode')

%% Projection

feature = 50;
digits = S * V'; % projection onto principal components: X = USV' --> U'X = SV'
TestMat = U(:,1:feature)'*test_images_reshape_demean; % PCA projection

%% build a linear classifier to identify two arbitrary digits

[threshold_arbitrary, w_arbitrary, sortdigitXs, sortdigitYs] = handmade_LDA(digits, labels, feature, 6, 8);
[threshold_hard, w_hard, sortdigitXs, sortdigitYs] = handmade_LDA(digits, labels, feature, 3, 5);
[threshold_easy, w_easy, sortdigitXs, sortdigitYs] = handmade_LDA(digits, labels, feature, 0, 1);

%% Quantify the accuracy on the test data

accuracy_handmade_LDA_arbitrary_training = accuracy_rate_handmade_LDA(w_arbitrary, digits(1:feature,:), labels, 6, 8, threshold_arbitrary);
accuracy_handmade_LDA_arbitrary = accuracy_rate_handmade_LDA(w_arbitrary, TestMat, test_labels, 6, 8, threshold_arbitrary);
accuracy_handmade_LDA_hard_training = accuracy_rate_handmade_LDA(w_hard, digits(1:feature,:), labels, 4, 9, threshold_hard);
accuracy_handmade_LDA_hard = accuracy_rate_handmade_LDA(w_hard, TestMat, test_labels, 4, 9, threshold_hard);
accuracy_handmade_LDA_easy_training = accuracy_rate_handmade_LDA(w_easy, digits(1:feature,:), labels, 0, 1, threshold_easy);
accuracy_handmade_LDA_easy = accuracy_rate_handmade_LDA(w_easy, TestMat, test_labels, 0, 1, threshold_easy);

%% Build a linear classifier to identify three arbitrarily picked digits

training_three_digits = digits(1:feature,:);
training_three_digits = digits(find(labels == 0 | labels == 1 | labels == 2));
training_three_digits_labels = labels(find(labels == 0 | labels == 1 | labels == 2));
test_three_digits = TestMat(find(test_labels == 0 | test_labels == 1 | test_labels == 2));
test_three_digits_labels = test_labels(find(test_labels == 0 | test_labels == 1 | test_labels == 2));
class = classify(test_three_digits, training_three_digits, training_three_digits_labels, 'linear');
num_of_correct = 0;
for i = 1:length(test_three_digits)
    if class(i) == test_three_digits_labels(i)
        num_of_correct = num_of_correct + 1;
    end
end
accuracy_three_digits_test = num_of_correct / length(test_three_digits);
class = classify(training_three_digits, training_three_digits, training_three_digits_labels, 'linear');
num_of_correct = 0;
for i = 1:length(test_three_digits)
    if class(i) == test_three_digits_labels(i)
        num_of_correct = num_of_correct + 1;
    end
end
accuracy_three_digits_training = num_of_correct / length(test_three_digits);

%% decision tree classifier

tree = fitctree(digits(1:feature,:)', labels);

%% predict using decision tree classifier

predict_labels_dtreeclassifier = predict(tree, TestMat');
dtreeclassifier_training_labels = predict(tree, digits(1:feature,:)');
performance_dtreeclassifier_training = classifier_performance(dtreeclassifier_training_labels, labels);
performance_dtreeclassifier_test = classifier_performance(predict_labels_dtreeclassifier, test_labels);
hardest_pair_of_digits_dtreeclassifier = accuracy_rate(predict_labels_dtreeclassifier, test_labels, 4, 9);
hardest_pair_of_digits_dtreeclassifier_training = accuracy_rate(dtreeclassifier_training_labels, labels, 4, 9);
easiest_pair_of_digits_dtreeclassifier = accuracy_rate(predict_labels_dtreeclassifier, test_labels, 0, 1);
easiest_pair_of_digits_dtreeclassifier_training = accuracy_rate(dtreeclassifier_training_labels, labels, 0, 1);

%% SVM classifier

Mdl = fitcecoc(digits(1:feature,:)', labels);

%% predict using SVM classifier

predict_labels = predict(Mdl, TestMat');
SVM_training_labels = predict(Mdl, digits(1:feature,:)');
performance_SVM_training = classifier_performance(SVM_training_labels, labels);
performance_SVM_test = classifier_performance(predict_labels, test_labels);
hardest_pair_of_digits = accuracy_rate(predict_labels, test_labels, 4, 9);
hardest_pair_of_digits_training = accuracy_rate(SVM_training_labels, labels, 4, 9);
easiest_pair_of_digits = accuracy_rate(predict_labels, test_labels, 0, 1);
easiest_pair_of_digits_training = accuracy_rate(SVM_training_labels, labels, 0, 1);

%% functions

function performance = classifier_performance(predict_digits, correct_digits)
    number_of_correct = 0;
    for j = 1:length(predict_digits)
        if predict_digits(j) == correct_digits(j)
            number_of_correct = number_of_correct + 1;
        end
    end
    performance = number_of_correct / length(predict_digits);
end

function accuracy_handmade_LDA = accuracy_rate_handmade_LDA(w, TestMat, test_labels, digitX, digitY, threshold)
    pval = w'*TestMat;
    test_digitXs = find(test_labels == digitX);
    test_digitYs = find(test_labels == digitY);
    numCorrect = 0;
    index = 0;
    for k = 1:length(test_digitXs)
        index = index + 1;
        if pval(test_digitXs(index)) < threshold;
            numCorrect = numCorrect + 1;
        end
    end
    index = 0;
    for k = 1:length(test_digitYs)
        index = index + 1;
        if pval(test_digitYs(index)) > threshold;
            numCorrect = numCorrect + 1;
        end
    end
    accuracy_handmade_LDA = numCorrect / (length(test_digitXs) + length(test_digitYs));
end
    

function [threshold,w,sortdigitXs,sortdigitYs] = handmade_LDA(digits, labels, feature, digitX, digitY)
    ndigitXs = length(find(labels == digitX));
    ndigitYs = length(find(labels == digitY));
    digitXs = digits(1:feature,find(labels == digitX));
    digitYs = digits(1:feature,find(labels == digitY));
    
    % Calculate scatter matrices

    mdigitXs = mean(digitXs,2);
    mdigitYs = mean(digitYs,2);

    Sw = 0; % within class variances
    for k = 1:ndigitXs
        Sw = Sw + (digitXs(:,k) - mdigitXs)*(digitXs(:,k) - mdigitXs)';
    end
    for k = 1:ndigitYs
        Sw =  Sw + (digitYs(:,k) - mdigitYs)*(digitYs(:,k) - mdigitYs)';
    end

    Sb = (mdigitXs-mdigitYs)*(mdigitXs-mdigitYs)'; % between class
    
    % Find the best projection line

    [V2, D] = eig(Sb,Sw); % linear disciminant analysis
    [lambda, ind] = max(abs(diag(D)));
    w = V2(:,ind);
    w = w/norm(w,2);
    
    % Project onto w

    vdigitXs = w'*digitXs;
    vdigitYs = w'*digitYs;
    
    % Make zeros below the threshold

    if mean(vdigitXs) > mean(vdigitYs)
        w = -w;
        vdigitXs = -vdigitXs;
        vdigitYs = -vdigitYs;
    end
    
    % Find the threshold value

    sortdigitXs = sort(vdigitXs);
    sortdigitYs = sort(vdigitYs);

    t1 = length(sortdigitXs);
    t2 = 1;
    while sortdigitXs(t1) > sortdigitYs(t2)
        t1 = t1 - 1;
        t2 = t2 + 1;
    end
    threshold = (sortdigitXs(t1) + sortdigitYs(t2))/2;
end

function accuracy = accuracy_rate(predict_labels, test_labels, digitX, digitY)
    predict_two_digits = predict_labels(find(test_labels == digitX | test_labels == digitY));
    test_two_digits = test_labels(find(test_labels == digitX | test_labels == digitY));
    nCorrect = 0;
    for i = 1:length(predict_two_digits)
        if predict_two_digits(i) == test_two_digits(i);
            nCorrect = nCorrect + 1;
        end
    end
    accuracy = nCorrect / length(predict_two_digits);
end

function images_reshape = reshapeimages(images)
    [m,~,num_of_images] = size(images);
    images_reshape = zeros(m^2, num_of_images);
    for k = 1:num_of_images
        images_reshape(:,k) = im2double(reshape(images(:,:,k),m^2,1));
    end
end