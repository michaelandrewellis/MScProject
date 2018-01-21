import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import normalize
from sklearn.metrics import classification_report, confusion_matrix, matthews_corrcoef, accuracy_score, make_scorer, \
    roc_curve, roc_auc_score
from sklearn.model_selection import cross_val_score, GridSearchCV, KFold
from scipy import interp
from sklearn.ensemble import VotingClassifier
from sklearn.utils import resample
from sklearn.svm import SVC
from sklearn.utils import shuffle
import pickle
from numpy import mean, std

from matplotlib import pyplot as plt
from matplotlib2tikz import save as tikz_save
from sklearn.model_selection import StratifiedKFold
from sklearn.feature_selection import RFECV

allele = 'fibroblasts'
positive_data = 'training_data_positive_binders_GRLCL.csv'
negative_data = 'training_data_negative_binders_50000_GRLCL.csv'


def cross_validation_with_sample_weights(n):
    # load data and set output values
    X_positive = pd.DataFrame.from_csv(positive_data)
    X_negative = pd.DataFrame.from_csv(negative_data)
    X_positive['output'] = 1
    X_negative['output'] = 0
    X_positive['sample_weight'] = 1
    X_negative['sample_weight'] = X_positive.__len__()/X_negative.__len__()


    X = pd.concat([X_positive, X_negative])
    X = shuffle(X)
    chunk_size = X.__len__()//n
    for i in range(0,10):
        test = X[chunk_size*i:chunk_size*(i+1)] # maybe fix this to include all examples
        train = X.drop(test.index)
        clf = RandomForestClassifier(n_estimators=1000, n_jobs=-1, max_features=20)
        
        

    y = X['output']
    X[['length1', 'distance']] = X[['length1', 'distance']].applymap(str)
    X.drop('output', axis=1, inplace=True)
    X = pd.get_dummies(X, dummy_na=True)

    # save columns to be able to recreate shape of data
    with open('columns_' + allele + '.txt', 'w') as f:
        for col in X.columns:
            f.write(col + '\n')
    MCC = make_scorer(matthews_corrcoef)
    clf = RandomForestClassifier(n_estimators=1000, n_jobs=-1, max_features=20)
    scores = cross_val_score(clf, X, y=y, cv=n, scoring=MCC)
    print(mean(scores))
    print(mean(scores) - 2 * std(scores), mean(scores) + 2 * std(scores))
    features = features_importance(clf, X, y)
    features = pd.DataFrame(features)
    features.columns = ['score', 'feature']
    new_features = pd.DataFrame()
    for i in range(1, 10):
        s = features[features.feature.str.startswith('seq_' + str(i))].score.sum()
        new_features = new_features.append([['seq_' + str(i), s]])
    for i in range(1, 5):
        for j in range(1, 8):
            s = features[features.feature.str.startswith('cleavage' + str(i) + '_' + str(j))].score.sum()
            new_features = new_features.append([['cleavage' + str(i) + '_' + str(j), s]])
    for i in range(1, 3):
        for j in range(1, 5):
            s = features[features.feature.str.startswith('binding' + str(i) + '_' + str(j))].score.sum()
            new_features = new_features.append([['binding' + str(i) + '_' + str(j), s]])
    # add feature importance for length1
    s = features[features.feature.str.startswith('length1')].score.sum()
    new_features = new_features.append([['length1', s]])
    # add feature importance for reversed
    s = features[features.feature.str.startswith('reversed')].score.sum()
    new_features = new_features.append([['reversed', s]])
    new_features.columns = ['feature', 'score']
    new_features.sort_values('score', ascending=False).to_csv('feature_scores.csv')

def cross_validation_final(n):
    # load data and set output values
    X_positive = pd.DataFrame.from_csv(positive_data)
    X_negative = pd.DataFrame.from_csv(negative_data)
    X_positive['output'] = 1
    X_negative['output'] = 0
    negative_weight = X_positive.__len__()/X_negative.__len__()

    # create balanced dataset
    X_negative = X_negative[0:X_positive.shape[0]]
    #X_negative = pd.concat([X_negative[0:55],X_negative[-414:]]) # For combined datasets?
    X = pd.concat([X_positive, X_negative])

    y = X['output']
    X[['length1', 'distance']] = X[['length1', 'distance']].applymap(str)
    X.drop('output', axis=1, inplace=True)
    X = pd.get_dummies(X, dummy_na=True)
    
    # save columns to be able to recreate shape of data
    with open('columns_' + allele + '.txt', 'w') as f:
        for col in X.columns:
            f.write(col + '\n')
    MCC = make_scorer(matthews_corrcoef)
    clf = RandomForestClassifier(n_estimators= 1000, n_jobs=-1, max_features= 20)#,class_weight={0:negative_weight,1:1})
    scores = cross_val_score(clf,X,y=y,cv=n,scoring=MCC)
    print(mean(scores))
    print(mean(scores)-2*std(scores),mean(scores)+2*std(scores))
    features = features_importance(clf,X,y)
    features = pd.DataFrame(features)
    features.columns = ['score','feature']
    features.sort('score').to_csv('feature_importances_dummies.csv')
    new_features = pd.DataFrame()
    for i in range(1,10):
        s = features[features.feature.str.startswith('seq_'+str(i))].score.sum()
        new_features = new_features.append([['seq_'+str(i),s]])
    for i in range(1,5):
        for j in range(1,8):
            s = features[features.feature.str.startswith('cleavage' + str(i)+'_'+str(j))].score.sum()
            new_features = new_features.append([['cleavage' + str(i)+'_'+str(j), s]])
    for i in range(1,3):
        for j in range(1,5):
            s = features[features.feature.str.startswith('binding' + str(i)+'_'+str(j))].score.sum()
            new_features = new_features.append([['binding' + str(i)+'_'+str(j), s]])
    # add feature importance for length1
    s = features[features.feature.str.startswith('length1')].score.sum()
    new_features = new_features.append([['length1', s]])
    # add feature importance for reversed
    s = features[features.feature.str.startswith('reversed')].score.sum()
    new_features = new_features.append([['reversed', s]])
    new_features.columns = ['feature','score']
    new_features.sort_values('score',ascending=False).to_csv('feature_scores.csv')
    
    

def setup_data():
    # load data and set output values
    X_positive = pd.DataFrame.from_csv(positive_data)
    X_negative = pd.DataFrame.from_csv(negative_data)
    X_positive['output'] = 1
    X_negative['output'] = 0

    # create balanced dataset
    X_negative = X_negative[0:X_positive.shape[0]]
    X = pd.concat([X_positive, X_negative])

    # split data into training, test and validation
    y = X['output']
    X[['length1', 'distance']] = X[['length1', 'distance']].applymap(str)
    X.drop('output', axis=1, inplace=True)
    if 'ic50' in X.columns:
        ic50 = X.ic50
        X = pd.get_dummies(X.drop('ic50', axis=1), dummy_na=True)
        X['ic50'] = ic50.map(lambda x: 1 - (np.log(x) / np.log(ic50.max())))
    else:
        X = pd.get_dummies(X, dummy_na=True)

    # save columns to be able to recreate shape of data
    with open('columns/' + allele + '.txt', 'w') as f:
        for col in X.columns:
            f.write(col + '\n')

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)
    return (X_train, y_train, X_test, y_test)


def train_test():
    X_train, y_train, X_test, y_test = setup_data()
    rf, rf_results, rf_score = tune_parameters(X_train, y_train, RandomForestClassifier(),
                                               {'max_features': [5, 10, 20, 50, 100, 200, None],
                                                'n_estimators': [1000],
                                                'n_jobs': [-1]})
    rf_results.to_csv('classifier_results/' + allele + '_rf.csv')
    nn, nn_results, nn_score = tune_parameters(X_train, y_train, MLPClassifier(),
                                               {'hidden_layer_sizes': [5, 10, 20, 50],
                                                'activation': ['identity', 'logistic', 'tanh', 'relu'],
                                                'max_iter': [1000]})
    nn_results.to_csv('classifier_results/' + allele + '_nn.csv')
    svm, svm_results, svm_score = tune_parameters(X_train, y_train, SVC(),
                                                  {'C': [0.1, 1, 10],
                                                   'kernel': ['linear', 'poly', 'rbf'],
                                                   'probability': [False]})
    svm_results.to_csv('classifier_results/' + allele + '_svm.csv')
    performance = []
    for model, name in zip([rf, nn, svm], ['Random Forest', 'Neural Network', 'SVM']):
        y_pred = model.predict(X_test)
        performance.append([name, matthews_corrcoef(y_test, y_pred)])
        cm = pd.DataFrame(confusion_matrix(y_test, y_pred))
        cm.columns = ['0', '1']
        cm.set_index(['0', '1'])
        cm.to_csv('classifier_results/' + allele + '_' + name + '_confusion_matrix.csv')
    pickle.dump(rf, open('classifiers/' + allele + '.pickle', 'wb'))
    pd.DataFrame(performance).to_csv('classifier_results/' + allele + '_test.csv')


def tune_parameters(X, y, model, params):
    clf = GridSearchCV(model, params, scoring='f1')
    clf.fit(X, y)
    rf_results = pd.DataFrame(clf.cv_results_)
    return (clf.best_estimator_, rf_results, clf.best_score_)


def cross_validation(X, y, model):
    # train classifiers
    scores = cross_val_score(model, X, y, cv=5, scoring='f1')
    return (scores)


def features_importance(model, X,y):
    model.fit(X,y)
    l = []
    for x, y in sorted(zip(model.feature_importances_, X.columns)):
        l.append([x, y])
    return(l)


def feature_selection():
    X, y, X_test, y_test = setup_data()
    model = RandomForestClassifier(n_estimators=500, n_jobs=-1)
    rfecv = RFECV(estimator=model, step=20, cv=StratifiedKFold(2),
                  scoring='f1')
    rfecv.fit(X, y)

    print("Optimal number of features : %d" % rfecv.n_features_)

    # Plot number of features VS. cross-validation scores
    for x, y in sorted(zip(rfecv.ranking_, X.columns)):
        print(x, y)
    plt.figure()
    plt.xlabel("Number of features selected")
    plt.ylabel("Cross validation score (f1)")
    plt.plot(range(1, len(rfecv.grid_scores_) * 20 + 1, 20), rfecv.grid_scores_)
    tikz_save('writeup/images/feature_selection_2.tex',
              figureheight='\\figureheight',
              figurewidth='\\figurewidth')
    plt.show()


def fix_columns(df, allele):
    cols = []
    with open('columns/' + allele + '.txt', 'r') as f:
        for line in f.readlines():
            cols.append(line.strip())
    X = pd.get_dummies(df, dummy_na=True)
    for col in X.columns:
        if col not in cols:
            X.drop(col, axis=1, inplace=True)
    for col in cols:
        if col not in X.columns:
            X[col] = 0
    X = X[cols]
    return (X)


def test_invented_binders(allele):
    df = pd.DataFrame.from_csv('training_data/test_data_negative_binders_' + allele + '.csv')
    df[['length1', 'distance']] = df[['length1', 'distance']].applymap(str)
    X = fix_columns(df, allele)
    clf = pickle.load(open('classifiers/' + allele + '.pickle', 'rb'))
    y_pred = clf.predict(X)
    y_true = [0] * len(y_pred)
    print(confusion_matrix(y_true, y_pred))


def test_general_classifier():
    df = pd.DataFrame.from_csv('training_data/training_data_with_cleavage_positive_C1R.csv')
    df[['length1', 'distance']] = df[['length1', 'distance']].applymap(str)
    X = fix_columns(df, 'non_C1R')
    clf = pickle.load(open('classifiers/non_C1R.pickle', 'rb'))
    y_pred = clf.predict(X)
    y_true = [0] * len(y_pred)
    print(accuracy_score(y_true, y_pred))


def plot_ROC_curve(n):
    # load data and set output values
    X_positive = pd.DataFrame.from_csv(positive_data)
    X_negative = pd.DataFrame.from_csv(negative_data)
    X_positive['output'] = 1
    X_negative['output'] = 0
    negative_weight = X_positive.__len__() / X_negative.__len__()
    
    scores=[]
    for i in range(n):
        # create balanced dataset
        #X_negative = resample(X_negative,replace=False,n_samples=X_positive.shape[0])#X_negative[0:X_positive.shape[0]]
        # X_negative = pd.concat([X_negative[0:55],X_negative[-414:]]) # For combined datasets?
        X = pd.concat([X_positive, X_negative])
    
        y = X['output']
        X[['length1', 'distance']] = X[['length1', 'distance']].applymap(str)
        X.drop('output', axis=1, inplace=True)
        X = pd.get_dummies(X, dummy_na=True)
        clf = RandomForestClassifier(n_estimators=1000, n_jobs=-1,
                                     max_features=20, class_weight={0:negative_weight,1:1})
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.2)
        clf.fit(X_train,y_train)
        y_score = clf.predict_proba(X_test)[:,1]
        roc_score = roc_auc_score(y_test,y_score)
        scores.append(roc_score)
    print(mean(scores),std(scores))
    fpr,tpr,_ = roc_curve(y_test,y_score)
    plt.plot(fpr, tpr, color='darkorange', label='ROC curve (area = %0.2f)' % roc_score)
    plt.show()

def plot_multiple_roc_curves():
    X_positive = pd.DataFrame.from_csv(positive_data)
    X_negative = pd.DataFrame.from_csv(negative_data)
    X_positive['output'] = 1
    X_negative['output'] = 0
    #X_negative = resample(X_negative,replace=False,n_samples=X_positive.shape[0])#X_negative[0:X_positive.shape[0]]
    X = pd.concat([X_positive, X_negative])
    y = X['output']
    X[['length1', 'distance']] = X[['length1', 'distance']].applymap(str)
    X.drop('output', axis=1, inplace=True)
    X = pd.get_dummies(X, dummy_na=True)
    
    kf = KFold(shuffle=True, n_splits=10)
    tprs = []
    auc_scores = []
    base_fpr = np.linspace(0, 1, 101)

    plt.figure(figsize=(5, 5))

    for i in kf.split(X):
        train, test = i[0],i[1]
        model = RandomForestClassifier(n_estimators=1000, n_jobs=-1,
                                     max_features=20)
        model.fit(X.iloc[train], y.iloc[train])
        y_score = model.predict_proba(X.iloc[test])
        fpr, tpr, _ = roc_curve(y.iloc[test], y_score[:, 1])
        plt.plot(fpr, tpr, 'b', alpha=0.15)
        tpr = interp(base_fpr, fpr, tpr)
        tpr[0] = 0.0
        tprs.append(tpr)
        auc_scores.append(roc_auc_score(y.iloc[test],y_score[:,1]))
    tprs = np.array(tprs)
    mean_tprs = tprs.mean(axis=0)
    std = tprs.std(axis=0)

    tprs_upper = np.minimum(mean_tprs + std, 1)
    tprs_lower = mean_tprs - std

    plt.plot(base_fpr, mean_tprs, 'b')
    plt.fill_between(base_fpr, tprs_lower, tprs_upper, color='grey', alpha=0.3)
    plt.show()
    print(mean(auc_scores))
    
if __name__ == '__main__':
    # tune_parameters('Random Forest')
    # feature_selection()
    # test_model()
    # oneclass()

    # test_invented_binders('C0702')
    #train_test()
    
    #plot_ROC_curve(10)
    plot_multiple_roc_curves()
    
    ##cross_validation_final(10)
    
    # test_general_classifier()
