import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

def calc_mse(data, fit):
    # 计算评估指标，这里可以使用均方误差、拟合优度或其他适当的指标
    mse = np.mean((data - stats.norm.pdf(data, *fit)) ** 2)
    return mse

def calc_mle(data, mean, sd):
    mle = np.sum(stats.norm.pdf(data, loc=mean, scale=sd))
    return mle


def simulate_data():
    # Generate sample data (replace with your own dataset)
    np.random.seed(0)
    data = np.concatenate([np.random.normal(0, 1, 500), 
                        np.random.normal(5, 1, 500)])
    return data

def measure_by_mle(data, parameter_file_name):
    data = np.array(data).reshape(-1, 1)

    # 初始化最优结果
    best_mle_score = float('-inf')
    best_fit_before_node = None
    best_fit_after_node = None
    best_node = None

    # 从左至右循环遍历所有节点位置
    # for node in np.unique(data):
    for node in np.linspace(np.min(data), np.max(data), 100):
        # 分段拟合
        data_before_node = data[data < node].reshape(-1, 1)
        data_after_node = data[data >= node].reshape(-1, 1)

        if len(data_before_node) <= 0:
            continue

        fit_before_node = stats.norm.fit(data_before_node)
        fit_after_node = stats.norm.fit(data_after_node)
        
        # try MLE
        # print("Likelihood", stats.norm.pdf(data_before_node, 
        #                             loc=fit_before_node[0], 
        #                             scale=fit_before_node[1]))
        # mle = np.sum(stats.norm.pdf(data_before_node, 
        #                             loc=fit_before_node[0], 
        #                             scale=fit_before_node[1])) + \
        #     np.sum(stats.norm.pdf(data_after_node, 
        #                             loc=fit_after_node[0], 
        #                             scale=fit_after_node[1]))
        mle = calc_mle(data_before_node, fit_before_node[0], fit_before_node[1]) + \
            calc_mle(data_after_node, fit_after_node[0], fit_after_node[1])

        #print(f"{node}\t{fit_before_node}\t{fit_after_node}\t{mse}")

        # 更新最优结果
        if mle > best_mle_score:
            best_mle_score = mle
            best_fit_before_node = fit_before_node
            best_fit_after_node = fit_after_node
            best_node = node

    # print("Best cutoff:", best_node)
    # print("Best mse:", best_mle_score)
    # print("Best Fit Before Node (Normal):", best_fit_before_node)
    # print("Best Fit After Node (Normal):", best_fit_after_node)

    #visual_optimal_cutoff(data, best_node, 
    #                      best_fit_before_node, best_fit_after_node)

    with open(parameter_file_name, 'w') as pf:
        # pf.write(str(optimal_cutoff))
        pf.write(f"Best cutoff: {best_node}\n")
        pf.write(f"Best mle: {best_mle_score}\n")
        pf.write(f"Best Fit Before Node (Normal): {best_fit_before_node}\n")
        pf.write(f"Best Fit After Node (Normal): {best_fit_after_node}")


def measure_by_mse(data):
    # 初始化最优结果
    best_score = float('inf')
    best_fit_before_node = None
    best_fit_after_node = None
    best_node = None

    # 从左至右循环遍历所有节点位置
    # for node in np.unique(data):
    for node in np.linspace(np.min(data), np.max(data), 100):
        # 分段拟合
        data_before_node = data[data < node].reshape(-1, 1)
        data_after_node = data[data >= node].reshape(-1, 1)

        if len(data_before_node) <= 0:
            continue

        fit_before_node = stats.norm.fit(data_before_node)
        fit_after_node = stats.norm.fit(data_after_node)

        # 计算评估指标，这里可以使用均方误差、拟合优度或其他适当的指标
        # mse = np.mean((data_before_node - 
        #                stats.norm.pdf(data_before_node, 
        #                               *fit_before_node)) ** 2) + \
        #       np.mean((data_after_node - 
        #                stats.norm.pdf(data_after_node, 
        #                               *fit_after_node)) ** 2)
        mse = calc_mse(data_before_node, fit_before_node) + \
            calc_mse(data_after_node, fit_after_node)

        #print(f"{node}\t{fit_before_node}\t{fit_after_node}\t{mse}")
        # 更新最优结果
        if mse < best_score:
            best_score = mse
            best_fit_before_node = fit_before_node
            best_fit_after_node = fit_after_node
            best_node = node


    print("Best cutoff:", best_node)
    print("Best mse:", best_score)
    print("Best Fit Before Node (Normal):", best_fit_before_node)
    print("Best Fit After Node (Normal):", best_fit_after_node)
    visual_optimal_cutoff(data, best_node, 
                          best_fit_before_node, best_fit_after_node)

def visual_optimal_cutoff(data, best_node, 
                          best_fit_before_node, best_fit_after_node):
    # 可视化最优拟合结果
    data_before_best_node = data[data < best_node].reshape(-1, 1)
    data_after_best_node = data[data >= best_node].reshape(-1, 1)

    plt.hist(data, bins=30, density=True, alpha=0.6, label="Data")
    x_before_best_node = np.linspace(np.min(data_before_best_node), 
                    np.max(data_before_best_node), 100)
    plt.plot(x_before_best_node, 
            stats.norm.pdf(x_before_best_node, *best_fit_before_node), 
            label="Before Node (Normal)")

    x_after_best_node = np.linspace(np.min(data_after_best_node), 
                    np.max(data_after_best_node), 100)
    plt.plot(x_after_best_node, 
            stats.norm.pdf(x_after_best_node, *best_fit_after_node), 
            label="After Node (Normal)")
    plt.legend()
    plt.show()

if __name__ == "__main__":
    data = simulate_data()
    measure_by_mle(data)