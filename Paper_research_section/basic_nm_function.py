#NM V.1.1
import matplotlib.pyplot as plt
import numpy as np
import sympy as sp
def graph_maker(func, plot=None, xl='x', yl='y', line=False):
    '''
    Function สำหรับวาดกราฟแบบขี้เกียจ รับค่าดังนี้
    
    func:             ฟังก์ชันที่ใช้
    plot              จุดที่จะ plot
    xl:               ชื่อแกน x
    yl:               ชื่อแกน y
    line:             สร้างเส้นตรง y = x
    '''
    x = np.linspace(plot[0], plot[1])
    y = func(x)
    plt.figure(figsize=(12,6));
    plt.grid();
    plt.xlabel(xl, fontsize=20);
    plt.ylabel(yl, fontsize=20);
    plt.plot(x, y, color='blue');
    if line == True:
        plt.plot(x, x, color='black');
    plt.show();

def root_finding_method(method, f, ini=None, e_s=0.0005, num=20):
    '''
    ฟังก์ชันสำหรับหาค่ารากแบบ NM ประกอบไปด้วยโหมดทั่งหมด 4 โหมดคือ
    method = 'Bisection Method'
    method = 'False-Position Method'
    method = 'Newton Raphson Method'
    method = 'Fixed-Point Iteration Method'
    และส่วนเสริมที่สำคัญดังนี้
    ini       ใช้สำหรับกำหนด x0 , xl และ xu
    ข้อควรระวัง Bisection Method ต้องเลือก f(xl) และ f(xu) ที่มีเครื่องหมายตรงกันข้ามกัน หรือเมื่อ f(xl)f(xu) < 0
    ข้อควรระวัง False-Position Method ต้องเลือก f(xl) และ f(xu) ที่มีเครื่องหมายตรงกันข้ามกัน หรือเมื่อ f(xl)f(xu) < 0
    ข้อควรระวัง Fixed-Point Iteration Method ต้องเริ่มจาก f(x) = 0 แล้วทำให้ x = g(x) ในการหา g(x) แล้วนำค่านี้ไปแทนใน f(x) 
    '''
    def func(x):
        return f(x) - x
    s = 50
    if method == 'Bisection Method':
        xl = ini[0]
        xu = ini[1]
        for i in range(num):
            xr = (xl + xu)/2
            if func(xr)*func(xl) < 0:
                xu = xr
            else:
                xl = xr
            if i > 1:
                epsilon = abs(func(xr))
                if epsilon < e_s:
                    print('-'*s)
                    print(f'Bisection Method Result: {xr:.10f}')
                    print('-'*s)
    elif method == 'False-Position Method':
        xl = ini[0]
        xu = ini[1]
        for i in range(num):
            xr = xu - (func(xu)*(xl - xu))/(func(xl) - func(xu))   
            if func(xr)*func(xl) < 0:
                xu = xr
            else:
                xl = xr
            if i > 1:
                epsilon = abs(func(xr))
                if epsilon < e_s:
                    print('-'*s)
                    print(f'False-Position Method Result: {xr:.10f}')
                    print('-'*s)
    elif method == 'Newton Raphson Method':
        x = sp.Symbol('x')
        df = sp.diff(func(x), x)
        funcp = sp.lambdify(x, df, 'numpy')
        for i in range(num):
            ini_next = ini - (func(ini)/funcp(ini))
            epsilon = abs(func(ini_next))
            if epsilon < e_s:
                print('-'*s)
                print(f'Newton Raphson Method Result: {ini_next:.10f}')
                print('-'*s)
            ini = ini_next
    elif method == 'Fixed-Point Iteration Method':
        for i in range(num):
            ini_next = func(ini)
            epsilon = abs(func(ini_next))
            if epsilon < e_s:
                print('-'*s)
                print(f'Fixed-Point Iteration Method Result: {ini_next:.10f}')
                print('-'*s)
            ini = ini_next
            
def gen_traj(func=None, seed=0, num=10):
    '''
    Generate trajectory
    Argumenets:
        func       A function in the iteration
        seed       Initial value in the iteration
        num        Number of iteration
    '''
    x = [None] * (num+1)
    x[0] = seed
    for k in range(num):
        x[k+1] = func(x[k])
    return x

def phaseline_plot(f, fp, sample=2, iternum=5, **kwargs):
    ext = (fp[-1] - fp[0])/(len(fp)-1)
    fpext = np.linspace(min(fp)-ext, max(fp)+ext, num=2)
    plt.figure(figsize=(12,6))
    plt.plot(fp, [0]*len(fp), marker='o', linewidth=2, color='black')
    plt.plot(fpext, [0]*2, color='black', linewidth=2)
    plt.xticks(fp)
    plt.yticks([])
    fp_w_lr = [min(fp)-ext] + fp + [max(fp)+ext]
    num_arrows = len(fp_w_lr)-1 
    arrows = []
    for k in range(num_arrows):
        a = list(np.linspace(fp_w_lr[k], fp_w_lr[k+1], num=sample+2))
        a.pop(0)
        a.pop(-1)
        arrows += a
    shapes = [None] * len(arrows)
    for k in range(len(arrows)):
        seed = arrows[k]
        term = gen_traj(f, seed=seed, num=iternum)[-1]
        shapes[k] = 1 if term - seed > 0 else -1   
    plt.quiver(arrows, [0]*len(arrows), shapes, [0]*len(arrows), **kwargs)
    plt.show()

def traj_plot(func, seed=None, seeds=None, num=5, yscale='linear', **kwargs):
    if seed is not None:
        seeds = [seed]
    t = list(range(num+1))
    plt.figure(figsize=(12,12))
    for s in seeds:
        trj = gen_traj(func, seed=s, num=num)
        plt.plot(t, trj, marker='o', label=f'seed = {s}', **kwargs)
    plt.xticks(t)
    plt.yscale(yscale)
    plt.grid()
    plt.legend()
    plt.show();