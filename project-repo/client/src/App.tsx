import React, { useState, ChangeEvent, FormEvent } from 'react';
import './index.css';

interface PredictResponse {
  prediction: string;
}

interface SimilaritySearchResponse {
  'similarity search results': [string, string, number][]; // [sequence, protein name, similarity score]
}

function App() {
  const [sequence, setSequence] = useState<string>('');
  const [prediction, setPrediction] = useState<string>('');
  const [similarityResults, setSimilarityResults] = useState<[string, string, number][]>([]); // Updated to include protein names
  const [expandedSequenceIndex, setExpandedSequenceIndex] = useState<number | null>(null);

  const handleInputChange = (event: ChangeEvent<HTMLInputElement>) => {
    setSequence(event.target.value);
  };

  const handleSubmit = async (event: FormEvent<HTMLFormElement>) => {
    event.preventDefault();

    try {
      const predictResponse = await fetch('http://127.0.0.1:5000/predict', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({ target_sequence: sequence }),
      });

      const predictData: PredictResponse = await predictResponse.json();
      setPrediction(predictData.prediction);

      const similarityResponse = await fetch('http://127.0.0.1:5000/similarity_search', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({ target_sequence: sequence }),
      });

      const similarityData: SimilaritySearchResponse = await similarityResponse.json();
      setSimilarityResults(similarityData['similarity search results']);

    } catch (error) {
      console.error('Error:', error);
    }
  };

  const handleSequenceClick = (index: number) => {
    setExpandedSequenceIndex(expandedSequenceIndex === index ? null : index);
  };

  return (
    <div className="flex flex-col items-center p-4 bg-gray-100 min-h-screen">
      <h1 className="text-3xl font-bold mb-4">Protein Sequence Input</h1>
      <form onSubmit={handleSubmit} className="bg-white p-6 rounded-lg shadow-lg w-full max-w-md">
        <label className="block text-lg font-medium mb-2">
          Protein Sequence:
          <input
            type="text"
            value={sequence}
            onChange={handleInputChange}
            required
            className="mt-1 block w-full p-2 border border-gray-300 rounded-md"
          />
        </label>
        <button
          type="submit"
          className="mt-4 px-4 py-2 bg-blue-500 text-white font-semibold rounded-lg hover:bg-blue-600"
        >
          Submit
        </button>
      </form>

      {prediction && (
        <div className="mt-6 bg-white p-4 rounded-lg shadow-md w-full max-w-md">
          <h2 className="text-2xl font-semibold mb-2">Prediction Result:</h2>
          <p>{prediction}</p>
        </div>
      )}

{similarityResults.length > 0 && (
  <div className="mt-6 bg-white p-4 rounded-lg shadow-md w-full max-w-md">
    <h2 className="text-2xl font-semibold mb-2">Similarity Search Results:</h2>
    <div className="overflow-x-auto">
      <ul className="list-disc pl-5">
        {similarityResults.map((result, index) => (
          <li key={index} className="mb-2">
            <div className="flex flex-col items-start space-y-2">
              <button
                className="text-blue-500 hover:underline"
                onClick={() => handleSequenceClick(index)}
              >
                {expandedSequenceIndex === index
                  ? result[0] // Show full sequence when clicked
                  : `${result[0].slice(0, 50)}...`} {/* Show truncated sequence by default */}
              </button>
              <p className="text-sm text-gray-600">{result[1]}</p> {/* Protein Name */}
              <p className="text-sm text-gray-600">{result[2].toFixed(2)}%</p> {/* Similarity Score */}
            </div>
          </li>
        ))}
      </ul>
    </div>
  </div>
)}

    </div>
  );
}

export default App;
