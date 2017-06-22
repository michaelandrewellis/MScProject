import pandas as pd
from sklearn.model_selection import train_test_split 
from sklearn.ensemble import RandomForestClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import classification_report,confusion_matrix, matthews_corrcoef
from sklearn.naive_bayes import BernoulliNB
from matplotlib import pyplot as plt

df1 = pd.read_csv('training_data_3_positive.csv',index_col=0)
df1 = df1[df1['length']==9].drop(['seq_10','seq_11','seq_12'],axis=1)
df1 = df1.dropna()
df1.drop_duplicates()
print(df1.shape)
df2 = pd.read_csv('training_data_3_negative.csv',index_col=0)
print(df2.shape)
df = pd.concat([df1,df2])
df['distance'] = df['distance']-1
print(df.shape)
df = df.drop(['reversed','length','c_term','n_term'],axis=1)  # distance was causing problems - 21-25 was too predictive
#df = df.sample(frac=1).reset_index(drop=True) # shuffle rows
y = df['output']
df = df.drop('output',axis=1)
#df = df.drop(['cleavage1_3','cleavage1_4','cleavage2_0','cleavage2_1','cleavage2_2','cleavage3_3','cleavage3_4','cleavage4_0','cleavage4_1','cleavage4_2'],axis=1)
#df_nan = df[['seq_10','seq_11','seq_12']]
#df_main = df.drop(['seq_10','seq_11','seq_12'],axis=1)
df_main = pd.get_dummies(df.applymap(str),dummy_na=False)
#print(df_main.shape)
#df_nan = pd.get_dummies(df_nan.applymap(str),dummy_na=False)
#X = pd.concat([df_main,df_nan],axis=1)
X = df_main

X_train, X_test, y_train, y_test = train_test_split(X, y)

mlp = RandomForestClassifier(n_estimators=100)
#mlp = MLPClassifier(hidden_layer_sizes=(13,13,13),max_iter=2000,activation='logistic')
#mlp = BernoulliNB()
mlp.fit(X_train,y_train)
predictions = mlp.predict(X_test)
prediction_probs = [x[0] for x in mlp.predict_proba(X_test)]
print(confusion_matrix(y_test,predictions))
print(classification_report(y_test,predictions))
print(matthews_corrcoef(y_pred=predictions,y_true=y_test))
correct = [1-abs(x[0]-x[1]) for x in zip(predictions, y_test)]
df = pd.DataFrame(sorted(zip(prediction_probs,predictions,y_test,correct)),columns=['prediction_probs','predictions','y','correct'])
df['decile'] = pd.qcut(df['prediction_probs'], 10, labels=False)
decile_means = df.groupby('decile')['correct'].mean()
plt.plot(decile_means)
plt.show()
#for i in sorted(zip([abs(sum(x)/float(len(x))) for x in mlp.coefs_[0]],list(X.columns))):
#    print(i)
#for x in sorted(zip(mlp.feature_importances_,list(X.columns))):
 #   print(x)