from fitter import Fitter
from scipy import stats
import matplotlib.pyplot as plt

# simulation of distribution
# reference: https://zhuanlan.zhihu.com/p/420047068
def simulate_distribution(data, plot_file_name, parameter_file_name):
    f = Fitter(data)
    f.fit()
    with open(parameter_file_name, 'w') as pf:
        pf.write("***************** fitting result: summary **************\n")
        pf.write(f.summary().to_string())
        pf.write("***************** fitting result: get_best **************\n")
        pf.write(str(f.get_best()))
    
    # Save distribution plot
    f.plot_pdf(Nbest=3, lw=1, method='sumsquare_error')
    plt.savefig(plot_file_name)

if __name__ == "__main__":
    # data generation
    data = stats.gamma.rvs(2, loc=1.5, scale=2, size=10000)
    simulate_distribution(data, "test.jpeg")