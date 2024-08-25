import pandas as pd
import numpy as np
from sklearn.preprocessing import OneHotEncoder
from sklearn.model_selection import train_test_split
from sklearn.svm import SVC
from sklearn.metrics import accuracy_score
from sklearn.metrics import  precision_score, recall_score, f1_score, roc_curve, auc, precision_recall_curve
from difflib import SequenceMatcher
import itertools


def calculate_ngram_frequencies(sequence, n):
    ngram_freq = {}
    for i in range(len(sequence) - n + 1):
        ngram = sequence[i:i+n]
        ngram_freq[ngram] = ngram_freq.get(ngram, 0) + 1
    return ngram_freq

def create_feature_vector(row, n, feature_vector):
    sequence = row['Sequence']
    ngram_freq = calculate_ngram_frequencies(sequence, n)
    # Initializing the feature vector with zeros
    vector = feature_vector.copy()
    # Updating the feature vector with the frequency of each n-gram
    for ngram, freq in ngram_freq.items():
        vector[ngram] = freq
    return pd.Series(vector)


def preprocess_data():
    df=pd.read_csv('./db.csv')
    df=df.iloc[:111]

    df=df.drop(['Entry','Reviewed', 'Entry Name', 'Gene Names', 'Gene Ontology (biological process)', 'Gene Ontology (molecular function)'], axis=1)

    for name in df['Protein names']:
        df['Protein names'] = np.where(df['Protein names'].str.contains('endolysin'), 'endolysin', df['Protein names'])
        df['Protein names'] = np.where(df['Protein names'].str.contains('Endolysin'), 'endolysin', df['Protein names'])
        df['Protein names'] = np.where(df['Protein names'].str.contains('Holin'), 'holin', df['Protein names'])
        df['Protein names'] = np.where(df['Protein names'].str.contains('holin'), 'holin', df['Protein names'])
    
    proteins=df
    sequences=df['Sequence']
    data = {'Sequence': df['Sequence'].tolist()}  # Extracting the 'Sequence' column as a list
    data_dict = {'Sequence': data['Sequence']}
    datF = pd.DataFrame(data_dict)
    n=3
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"  # A string representing the 20 standard amino acids
    feature_vector = {"".join(aa): 0 for aa in itertools.product(amino_acids, repeat=n)}
    feature_vectors = datF.apply(lambda row: create_feature_vector(row, n, feature_vector), axis=1)
    result_df = pd.concat([datF['Sequence'], feature_vectors], axis=1)
    df = pd.concat([df, result_df.drop(columns=['Sequence'])], axis=1, join='inner')
    df=df.drop(['Sequence', 'Organism'], axis=1)
    enc=OneHotEncoder(sparse_output=False, handle_unknown='ignore')
    df['Protein names']=enc.fit_transform(np.array(df['Protein names']).reshape(-1,1))
    return df, proteins, sequences, feature_vector



def train_model(df):
    x=df.drop('Protein names', axis=1)
    y=df['Protein names']
    X_train, X_test, y_train, y_test = train_test_split(x, y, test_size=0.4, random_state=42)
    svc=SVC()
    svc.fit(X_train, y_train)
    y_pred=svc.predict(X_test)
    return y_pred, y_test, svc


def predict(prediction_sequence, seq_length, svc):
    ngram_freq = calculate_ngram_frequencies(prediction_sequence, 3)
    vector = feature_vector.copy()
    # Updating the feature vector with the frequency of each n-gram
    for ngram, freq in ngram_freq.items():
        vector[ngram] = freq
    pred_vec= pd.Series(vector)
    pred_vec=pd.DataFrame(np.array(pred_vec).reshape(1,-1))
    length=pd.DataFrame(np.array([seq_length]))
    pred_vec=pd.concat([pred_vec, length], axis="columns")
    result=svc.predict(pred_vec)
    if result[0]==1:
        protein="endolysin"
    else:
        protein="holin"
    return protein



def similarity_search(target_sequence, protein):
    a=proteins.loc[proteins['Protein names']==protein]
    sequences_list=a['Sequence']
    sequences_list=sequences_list.to_numpy()
    def similarity(seq1, seq2):
        """Calculate similarity score using SequenceMatcher."""
        return SequenceMatcher(None, seq1, seq2).ratio()

    # Calculate similarity for each sequence in the list
    similarities = [(seq, similarity(target_sequence, seq) * 100) for seq in sequences_list]
    
    # Sort the list based on similarity scores in descending order
    similarities.sort(key=lambda x: x[1], reverse=True)
    
    # Return the top N most similar sequences with their percentage match
    return similarities[:5]


df, proteins, sequences, feature_vector=preprocess_data()
y_pred, y_test, svc=train_model(df)

prediction_sequence="MKTKKQALKWILNTIGQGIDWDKMYGFQCMDLVVAYLYYVTDGKIAMWGNAIDAPKNNFKGTAKVIKNYPAFRPEEGDIVVWSYGNFSTYGHIAVVIDGDPYGDLQYITVAEQNWNGLGLYKQEVTTKRIHNYDGVSHFIRPKFKKTAKKEDNTPTKEKNNKKTKGKKLKVSTQRINYTMDKRGYKPKFVVIHNDAGSSSAQQYEQGLKNAGYSRYAQGVAHAYASDGYVWEAISEDRIAWHTGDGTNPGTGNFEGYGIEVCQSLGDRNTFLKNEQTVFQFIAEKLQKWNLPANRNTIRLHNEFIQTECPHASAYYHAGMNTKVDAYTKERQLKIKDYFIKQIRAYMKGSTPKSTVVKSSKSSGSLPKKKGQTSKSNIGKTFDFNGLSINVWGTKWYYENNTFTCNARQGIITRVGSPFTTAPQAGVLFYGQTVTYNQVAVNPKEPFVWISWITNNGTEVWMPIEVLDSNNKIIEQWGTFGW"
seq_length=len(prediction_sequence)
protein = predict("", seq_length, svc)

result = similarity_search(prediction_sequence, protein)






