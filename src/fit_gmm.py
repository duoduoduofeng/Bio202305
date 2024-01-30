import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture
import os

def fit_gmm(ori_data, plot_file_name, parameter_file_name):
    species_name = os.path.basename(plot_file_name).split(".")[0]

    # Create a GMM object
    n_components = 2  # Number of components/clusters
    gmm = GaussianMixture(n_components=n_components)

    # Fit the GMM to the data
    data = np.array(ori_data)
    gmm.fit(data.reshape(-1, 1))

    # Set up the figure and axis
    plt.clf()
    fig, ax = plt.subplots()

    # Generate x values to evaluate the GMM
    x = np.linspace(60, 100, 500)

    # Plot the GMM curve
    y = np.exp(gmm.score_samples(x.reshape(-1, 1)))
    ax.plot(x, y, label='GMM Curve')

    # Plot the data histogram
    ax.hist(data, bins=30, density=True, alpha=0.5, label='Data Histogram')
    ax.set_title(f'Gaussian Mixture Model for {species_name}')
    ax.set_xlabel('Similarity')
    ax.set_ylabel('Density')
    ax.legend()
    plt.savefig(plot_file_name)

    # Organize the parameters into a dictionary
    parameters = {
        "weights": gmm.weights_,
        "means": gmm.means_, 
        "covariances": gmm.covariances_
    }
    # Save the parameters to a file
    with open(parameter_file_name, 'w') as pf:
        pf.write(str(parameters))
        pf.write("\n")
    # np.save(parameter_file_name, parameters)

if __name__ == "__main__":
    # data generation
    data = stats.gamma.rvs(2, loc=1.5, scale=2, size=10000)
    fit_gmm(data, "test_gmm.jpeg", "test_gmm.parameters")