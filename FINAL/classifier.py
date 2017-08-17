import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split 
from sklearn.ensemble import RandomForestClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import normalize
from sklearn.metrics import classification_report,confusion_matrix, matthews_corrcoef,accuracy_score
from sklearn.model_selection import cross_val_score, GridSearchCV
from sklearn.ensemble import VotingClassifier
from sklearn.svm import SVC
import pickle

from matplotlib import pyplot as plt
from sklearn.model_selection import StratifiedKFold
from sklearn.feature_selection import RFECV

allele = 'GRLCL'
positive_data = 'training_data/training_data_with_cleavage_positive_GRLCL.csv'
negative_data = 'training_data/training_data_with_cleavage_negative.csv'

def setup_data():
    # load data and set output values
    X_positive = pd.DataFrame.from_csv(positive_data)
    X_negative = pd.DataFrame.from_csv(negative_data)
    X_positive['output'] = 1
    X_negative['output'] = 0
    
    # create balanced dataset
    X_negative = X_negative[0:X_positive.shape[0]]
    X = pd.concat([X_positive,X_negative])
    
    # split data into training, test and validation
    y = X['output']
    X[['length1','distance']] = X[['length1','distance']].applymap(str)
    X.drop('output',axis=1,inplace=True)
    X = pd.get_dummies(X, dummy_na=True)
    
    # save columns to be able to recreate shape of data
    with open('columns/'+allele+'.txt','w') as f:
        for col in X.columns:
            f.write(col+'\n')

    X, X_test, y, y_test = train_test_split(X, y,test_size=0.2)
    #X_train, X_validation, y_train, y_validation = train_test_split(X, y,test_size=0.2)
    return(X,y,X_test,y_test)

def train_test():
    X_train, y_train,X_test, y_test = setup_data()
    rf,rf_results,rf_score = tune_parameters(X_train,y_train,RandomForestClassifier(),
                                         {'max_features':[5,10,20,50,100,200,None],
                                          'n_estimators': [1000], 
                                          'n_jobs':[-1]})
    nn,nn_results,nn_score = tune_parameters(X_train,y_train,MLPClassifier(),
                                         {'hidden_layer_sizes':[5,10,20,50],
                                          'activation':['identity', 'logistic', 'tanh', 'relu'],
                                          'max_iter':[1000]})
    svm,svm_results,svm_score = tune_parameters(X_train,y_train,SVC(), 
                                           {'C':[0.1,1,10],
                                            'kernel': ['linear', 'poly', 'rbf'],
                                            'probability':[False]})
    rf_results.to_csv('classifier_results/'+allele+'_rf.csv')
    nn_results.to_csv('classifier_results/'+allele+'_nn.csv')
    svm_results.to_csv('classifier_results/'+allele+'_svm.csv')
    performance = []
    for model,name in zip([rf,nn,svm],['Random Forest','Neural Network','SVM']):
        y_pred = model.predict(X_test)
        performance.append([name,matthews_corrcoef(y_test,y_pred)])
        cm = pd.DataFrame(confusion_matrix(y_test,y_pred))
        cm.columns = ['0','1']
        cm.set_index(['0','1'])
        cm.to_csv('classifier_results/'+allele+'_'+name+'_confusion_matrix.csv')
    pickle.dump(rf,open('classifiers/'+allele+'.pickle','wb'))
    pd.DataFrame(performance).to_csv('classifier_results/'+allele+'_test.csv')
    
def tune_parameters(X,y,model,params):
    clf = GridSearchCV(model,params)
    clf.fit(X,y)
    rf_results = pd.DataFrame(clf.cv_results_)
    return(clf.best_estimator_,rf_results,clf.best_score_)

def cross_validation(X,y,model):
    # train classifiers
    scores = cross_val_score(model, X, y, cv=5, scoring='f1')
    return(scores)

def features_importance(model,X):
    model.fit(X)
    for x, y in sorted(zip(model.feature_importances_, X.columns)):
        print(x, y)   
        
def feature_selection():
    X,y,X_test,y_test = setup_data()
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
    plt.ylabel("Cross validation score (nb of correct classifications)")
    plt.plot(range(1, len(rfecv.grid_scores_)*20 + 1,20), rfecv.grid_scores_)
    plt.show()

def fix_columns(df,allele):
    cols = []
    with open('columns/'+allele+'.txt','r') as f:
        for line in f.readlines():
            cols.append(line.strip())
    X = pd.get_dummies(df,dummy_na=True)
    for col in X.columns:
        if col not in cols:
            X.drop(col,axis=1,inplace=True)
    for col in cols:
        if col not in X.columns:
            X[col]=0
    X=X[cols]
    return(X)
    
    

def test_invented_binders(allele):
    df= pd.DataFrame.from_csv('training_data/test_data_negative_binders_'+allele+'.csv')
    df[['length1', 'distance']] = df[['length1', 'distance']].applymap(str)
    X = fix_columns(df,allele)
    clf = pickle.load(open('classifiers/'+allele+'.pickle','rb'))
    y_pred = clf.predict(X)
    y_true = [0]*len(y_pred)
    print(accuracy_score(y_true,y_pred))

def test_general_classifier():
    df = pd.DataFrame.from_csv('training_data/training_data_with_cleavage_positive_C1R.csv')
    df[['length1', 'distance']] = df[['length1', 'distance']].applymap(str)
    X = fix_columns(df, 'non_C1R')
    clf = pickle.load(open('classifiers/non_C1R.pickle', 'rb'))
    y_pred = clf.predict(X)
    y_true = [0] * len(y_pred)
    print(accuracy_score(y_true, y_pred))
    
    
if __name__ == '__main__':
    
    #tune_parameters('Random Forest')
    #feature_selection()
    #test_model()
    #oneclass()
    
    #test_invented_binders('C0702')
    train_test()
    #test_general_classifier()
