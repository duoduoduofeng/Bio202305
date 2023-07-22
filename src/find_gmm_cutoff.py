import numpy as np
from sklearn.mixture import GaussianMixture

def find_gmm_cutoff(ori_data, parameter_file_name):
    # Reshape data to a column vector if necessary
    data = np.array(ori_data).reshape(-1, 1)
    
    # Fit a Gaussian Mixture Model with two components
    n_components = 2
    gmm = GaussianMixture(n_components=n_components)
    gmm.fit(data)

    # Calculate the posterior probabilities
    posterior_probs = gmm.predict_proba(data)
    #print(posterior_probs)

    # Find the optimal cutoff by maximizing the likelihood
    max_likelihood = -np.inf
    optimal_cutoff = None

    # The first two set were even more wrong.
    # for cutoff in np.linspace(np.min(data), np.max(data), 100):
    # for cutoff in np.unique(data):
    for cutoff in data:
        # Assign data points to components based on the posterior probabilities
        #left_component = data[posterior_probs[:, 0] > posterior_probs[:, 1]]
        #print(left_component)
        #right_component = data[posterior_probs[:, 0] <= posterior_probs[:, 1]]
        
        left_component = data[data < cutoff].reshape(-1, 1)
        right_component = data[data >= cutoff].reshape(-1, 1)
        #print(right_component)

        # filter the first data point
        if not len(left_component):
            continue

        # Calculate the log-likelihood of each component
        left_likelihood = gmm.score_samples(left_component).sum()
        right_likelihood = gmm.score_samples(right_component).sum()

        # Calculate the total likelihood
        total_likelihood = left_likelihood + right_likelihood

        # Update the maximum likelihood and cutoff if necessary
        if total_likelihood > max_likelihood:
            max_likelihood = total_likelihood
            optimal_cutoff = cutoff

    # print("Optimal Cutoff:", optimal_cutoff)
    
    # Save the parameters to a file
    with open(parameter_file_name, 'w') as pf:
        pf.write(str(optimal_cutoff))


if __name__ == "__main__":
    # Generate sample data (replace with your own dataset)
    np.random.seed(0)
    data = np.concatenate([np.random.normal(0, 1, 500), np.random.normal(5, 1, 500)])

    find_gmm_cutoff(data, "gmm_cutoff_sample.parameters")