import pandas as pd
from math import sqrt

def test_MCC_table(fout):
    datasets  = ['A0101','A0301','B0702','B2705','C0702','GRLCL','C1R','fibroblasts']
    df=pd.DataFrame()
    test_sizes=[]
    for dataset in datasets:
        test_sizes.append(pd.DataFrame.from_csv('classifier_results/'+dataset+'_Random Forest_confusion_matrix.csv').sum().sum())
        df = pd.concat([df,pd.DataFrame.from_csv('classifier_results/'+dataset+'_test.csv').set_index('0').round(3)],axis=1)
    df.columns = datasets
    df.index.name = None
    df = pd.concat([pd.DataFrame([test_sizes],columns=datasets),df])
    df=df.T
    df = df.rename(columns={0:'N'})
    df['N'] = df['N'].astype(int)
    df.index.name = 'Dataset'
    df.reset_index(inplace=True)
    df.columns = ['Dataset','N','RF','NN','SVM']
    with open(fout,'w') as f:
        f.write(df.to_latex(index=False,column_format='lcccc'))
        
def CM_table():
    pass

def CV_table(dataset,classifier,fout):
    df = pd.DataFrame.from_csv('classifier_results/'+dataset+'_'+classifier+'.csv')
    cols = [col for col in df.columns if col.startswith('param_')] + ['mean_test_score','std_test_score']
    df = df[cols]
    if classifier == 'rf':
        df = df.drop(['param_n_estimators','param_n_jobs'],axis=1)
    if classifier == 'svm':
        df = df.drop(['param_probability'],axis=1)
    if classifier == 'nn':
        df = df.drop(['param_max_iter'],axis=1)
    df.columns = [' '.join([x.capitalize() for x in col.split('_')]) for col in df.columns]
    df.columns = [col.replace('Param ', '') for col in df.columns]
    with open(fout,'w') as f:
        f.write(df.to_latex(index=False,column_format=len(df.columns)*'c'))
        
def all_CV_tables():
    for dataset in ['A0101','A0301','B0702','B2705','C0702','GRLCL','C1R','fibroblasts']:
        for classifier in ['rf','nn','svm']:
            CV_table(dataset,classifier,'writeup/tables/CV_table_'+dataset+'_'+classifier+'.tex')
            
def sensitivity_specificity(fout):
    datasets = ['GRLCL','B0702', 'C0702','C1R']
    table=[]
    for set in datasets:
        df = pd.DataFrame.from_csv('classifier_results/'+set+'_Random Forest_confusion_matrix.csv')
        tn = df.iloc[0,0]
        fn = df.iloc[1,0]
        fp = df.iloc[0,1]
        tp = df.iloc[1,1]
        sp = tn/(fp+tn)
        sn = tp/(tp+fn)
        pr = tp/(tp+fp)
        n = tp+fp+tn+fn
        acc = (tp+tn)/(tp+fp+tn+fn)
        mcc = (tp*tn-fp*fn)/(sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
        table.append([n,sp,sn,pr,acc,mcc])
    df = pd.DataFrame(table,columns = ['N','Specificity','Sensitivity','Precision','Accuracy','MCC'],index=datasets).round(3)
    df.index.name = 'Dataset'
    df.reset_index(inplace=True)
    df.N = df.N.astype(int)
    with open(fout,'w') as f:
        f.write(df[['Dataset','N','Sensitivity','Specificity','MCC']].to_latex(index=False,column_format='lcccc'))
        
    
    
        
if __name__ == '__main__':
    test_MCC_table('writeup/tables/test_MCC.tex')
    #all_CV_tables()
    #sensitivity_specificity('writeup/tables/sensitivity_specificity.tex')
    