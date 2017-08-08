import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split 
from sklearn.ensemble import RandomForestClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import normalize
from sklearn.metrics import classification_report,confusion_matrix, matthews_corrcoef
from sklearn.naive_bayes import BernoulliNB
from sklearn.model_selection import cross_val_score
from sklearn import svm
from matplotlib import pyplot as plt

df1 = pd.read_csv('training_data_3_positive_with_binding.csv',index_col=0)
df1 = df1[df1['length']==9].drop(['seq_10','seq_11','seq_12'],axis=1)
df1 = df1.dropna()
df1 = df1.drop_duplicates()
print(df1.shape)
df2 = pd.read_csv('training_data_3_negative_with_binding.csv',index_col=0)
print(df2.shape)
df = pd.concat([df1,df2])
print(df.columns)
df['distance'] = df['distance']-1
print(df.shape)
df = df.drop(['length2','length','c_term','n_term','distance'],axis=1)  # distance was causing problems - 21-25 was too predictive
#df = df.sample(frac=1).reset_index(drop=True) # shuffle rows
y = df['output']
df = df.drop('output',axis=1)
'''df = df.drop(['cleavage1_0', 'cleavage1_1',
       'cleavage1_3', 'cleavage1_4', 'cleavage2_0', 'cleavage2_1',
       'cleavage2_2', 'cleavage2_4', 'cleavage3_0',
       'cleavage3_1', 'cleavage3_3', 'cleavage3_4',
       'cleavage4_0', 'cleavage4_1', 'cleavage4_2',
       'cleavage4_4'],axis=1)'''
#df = pd.concat([df,df.bind1+df.bind2,df.seq_2+df.seq_9,df.bind1+df.cleavage3_2,df.bind1+df.cleavage2_4,df.seq_1+df.seq_2,df.cleavage3_2+df.bind2,df.cleavage1_2+df.seq_1,df.seq_9+df.cleavage4_3],axis=1)
#df = df.drop(['cleavage1_3', 'cleavage1_4','cleavage4_0', 'cleavage4_1','cleavage4_2'],axis=1)
#df = df.drop(['seq_3','seq_4','seq_5','seq_6','seq_7','seq_8'],axis=1)
#df_nan = df[['seq_10','seq_11','seq_12']]
#df_main = df.drop(['seq_10','seq_11','seq_12'],axis=1)
df_main = pd.concat([pd.get_dummies(df.drop('ic50',axis=1),dummy_na=False),df.ic50.map(lambda x: 1-np.log(x)/np.log(50000))#/df.ic50.map(np.log).max()
],axis=1) # CHOOSE IC50 TRANSFORMATION
#print(df_main.shape)
#df_nan = pd.get_dummies(df_nan.applymap(str),dummy_na=False)
#X = pd.concat([df_main,df_nan],axis=1)
X = df_main

X_train, X_test, y_train, y_test = train_test_split(X, y,test_size=0.2)
'''
X, X_test, y, y_test = train_test_split(X, y,test_size=0.2)
X_train, X_validation, y_train, y_validation = train_test_split(X, y,test_size=0.2)
'''

#mlp = RandomForestClassifier(n_estimators=500)
#mlp = MLPClassifier(hidden_layer_sizes=(10),activation='logistic')
#mlp = BernoulliNB()
mlp = svm.SVC(kernel='linear',probability=True)
scores = cross_val_score(mlp, X, y, cv=5,scoring='f1')
print("Accuracy: %0.2f (+/- %0.2f)" % (scores.mean(), scores.std() * 2))
mlp.fit(X_train,y_train)
predictions = mlp.predict(X_test)
print(confusion_matrix(y_test,predictions))
print(classification_report(y_test,predictions))
print(matthews_corrcoef(y_pred=predictions,y_true=y_test))
#correct = [1-abs(x[0]-x[1]) for x in zip(predictions, y_test)]
#prediction_probs = [x[0] for x in mlp.predict_proba(X_test)]
#df = pd.DataFrame(sorted(zip(prediction_probs,predictions,y_test,correct)),columns=['prediction_probs','predictions','y','correct'])
#df['decile'] = pd.qcut(df['prediction_probs'], 10, labels=False)
#decile_means = df.groupby('decile')['correct'].mean()
#plt.plot(decile_means)
#plt.show()
#for i in sorted(zip([abs(sum(x)/float(len(x))) for x in mlp.coefs_[0]],list(X.columns))):
 #   print(i)
#for x in sorted(zip(mlp.feature_importances_,list(X.columns))):
 #  print(x)
for coef,feature in sorted(zip(mlp.coef_[0],X.columns)):
    print(coef,feature)