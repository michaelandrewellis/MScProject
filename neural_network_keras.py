from keras.models import Sequential
from keras.layers import Dense
from sklearn.model_selection import train_test_split
from sklearn.metrics import precision_recall_fscore_support
from sklearn.utils import shuffle
import pandas as pd
import numpy
# fix random seed for reproducibility
numpy.random.seed(9)

# get data
df1 = pd.read_csv('training_data_3_positive.csv',index_col=0)
df1 = df1[df1['length']==9].drop(['seq_10','seq_11','seq_12'],axis=1)
df1 = df1.dropna()
df1.drop_duplicates()
df2 = pd.read_csv('training_data_3_negative.csv',index_col=0)
df2 = shuffle(df2,n_samples=df1.shape[0])
df = pd.concat([df1,df2])
df = df.drop(['reversed','length','distance','c_term','n_term'],axis=1)  # distance was causing problems - 21-25 was too predictive
y = df['output']
df = df.drop('output',axis=1)
#df = df.drop(['cleavage1_3','cleavage1_4','cleavage2_0','cleavage2_1','cleavage2_2','cleavage3_3','cleavage3_4','cleavage4_0','cleavage4_1','cleavage4_2'],axis=1)
#df_nan = df[['seq_10','seq_11','seq_12']]
#df_main = df.drop(['seq_10','seq_11','seq_12'],axis=1)
df_main = pd.get_dummies(df.applymap(str),dummy_na=False)
#df_nan = pd.get_dummies(df_nan.applymap(str),dummy_na=False)
#X = pd.concat([df_main,df_nan],axis=1)
X = df_main

# split data
X_train, X_test, y_train, y_test = train_test_split(X, y)

# create model
model = Sequential()
model.add(Dense(2, input_dim=636, activation='relu'))
#model.add(Dense(5, activation='relu'))
model.add(Dense(1, activation='sigmoid'))

# Compile model
model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['accuracy'])

# Fit the model
model.fit(X_train.as_matrix(), y_train.as_matrix(), epochs=6, batch_size=32,verbose=1)

# evaluate the model
scores = model.evaluate(X_test.as_matrix(), y_test.as_matrix())
y_pred = [round(x[0]) for x in model.predict(X_test.as_matrix()).tolist()]
y_prob = pd.DataFrame(model.predict_proba(X_test.as_matrix()))
y_prob.columns = ['pred']
print(y_prob.sort_values(by='pred'))
print(y_test)
print(precision_recall_fscore_support(y_test.tolist(),y_pred))
print("\n%s: %.2f%%" % (model.metrics_names[1], scores[1]*100))