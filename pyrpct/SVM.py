# import packages
from sklearn import svm
from joblib import dump, load
from Read import read_feature
from sklearn.model_selection import cross_val_predict
from sklearn.model_selection import GridSearchCV


# train
def svm_train(file, c_number, gamma, out):
    train_data, train_label = read_feature(file)
    my_svm = svm.SVC(kernel='rbf', C=c_number, gamma=gamma)
    model = my_svm.fit(train_data, train_label)
    dump(model, out)
    return train_label, train_data


# predict
def svm_predict(file, model, cv):
    test_data, test_label = read_feature(file)
    model = load(model)
    predict_label = cross_val_predict(model, test_data, test_label, cv=cv)
    return test_label, predict_label


# evaluate
def svm_evaluate(file, c_number, gamma, cv):
    test_data, test_label = read_feature(file)
    my_svm = svm.SVC(kernel='rbf', C=c_number, gamma=gamma)
    predict_label = cross_val_predict(my_svm, test_data, test_label, cv=cv)
    return test_label, predict_label


# grid search
def svm_grid(file):
    my_svm = svm.SVC(decision_function_shape="ovo", random_state=0)
    c_number = []
    for i in range(-5, 15 + 1, 2):
        c_number.append(2 ** i)
    gamma = []
    for i in range(-15, 3 + 1, 2):
        gamma.append(2 ** i)
    parameters = {'C': c_number, 'gamma': gamma}
    new_svm = GridSearchCV(my_svm, parameters, cv=5, scoring="accuracy", return_train_score=False, n_jobs=1)
    train_data, train_label = read_feature(file)
    model = new_svm.fit(train_data, train_label)
    best_c = model.best_params_['C']
    best_g = model.best_params_['gamma']
    return best_c, best_g


# acc
def svm_acc(y_t, y_p):
    tp = 0
    tn = 0
    fp = 0
    fn = 0
    y_t = list(y_t)
    y_p = list(y_p)
    for i in range(len(y_t)):
        t_label = y_t[i]
        p_label = y_p[i]
        if t_label == p_label == 0:
            tp += 1
        if t_label == p_label == 1:
            tn += 1
        if t_label == 0 and p_label == 1:
            fn += 1
        if t_label == 1 and p_label == 0:
            fp += 1
    acc = (tp + tn) / (tp + tn + fp + fn)
    print('\nPredict result: ' + str(round(acc, 3)*100) + '%')
