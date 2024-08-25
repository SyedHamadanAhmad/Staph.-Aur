# Filename - server.py

# Import flask and datetime module for showing date and time
from flask import Flask, request, jsonify
from model_final import preprocess_data, train_model, predict, similarity_search
from flask_cors import CORS



# Initializing flask app
app = Flask(__name__)
CORS(app)

# Preprocess data and train the model
df, proteins, sequences, feature_vector = preprocess_data()
y_pred, y_test, svc = train_model(df)


cache = {}

@app.route('/predict', methods=['POST'])
def predict_protein():
    data = request.json
    target_sequence = data.get('target_sequence', '')
    
    if not target_sequence:
        return jsonify({'error': 'No target sequence provided'}), 400
    
    # Check if the result is already in the cache
    else:
        # If not cached, compute the protein value
        protein = predict(target_sequence, len(target_sequence), svc)
        cache[target_sequence] = protein  # Store the result in the cache
        
        
    
    return jsonify({'prediction': protein})


@app.route('/similarity_search', methods=['POST'])
def similarity_serach():
    data=request.json
    target_sequence=data.get('target_sequence', '')
    protein=cache[target_sequence]
    similar_proteins=similarity_search(target_sequence, protein)
    return jsonify({'similarity search results': similar_proteins})


# Running app
if __name__ == '__main__':
    app.run(debug=True)
