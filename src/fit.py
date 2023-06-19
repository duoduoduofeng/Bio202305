from fitter import Fitter
from scipy import stats

# simulation of distribution
# reference: https://zhuanlan.zhihu.com/p/420047068
def simulate_distribution(data):
    f = Fitter(data)
    f.fit()
    print("***************** fitting result: summary **************")
    print(f.summary())
    print("***************** fitting result: get_best **************")
    print(f.get_best())
    print("***************** fitting result: plot_pdf **************")
    f.plot_pdf(Nbest=5, lw=2, method='sumsquare_error')
    print("***************** fitting result: fd_errors **************")
    f.df_errors()


if __name__ == "__main__":
    # data generation
    data = stats.gamma.rvs(2, loc=1.5, scale=2, size=10000)
    simulate_distribution(data)