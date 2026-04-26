import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

x = np.linspace(-2, 2, 1000)
y_cos = np.cos(x)
x_circle = np.linspace(-1, 1, 500)
y_circle_upper = 1 + np.sqrt(1 - x_circle**2)  # верхняя полуокружность
y_circle_lower = 1 - np.sqrt(1 - x_circle**2)  # нижняя полуокружность

ax1 = axes[0]
ax1.plot(x, y_cos, 'b-', linewidth=2)
ax1.plot(x_circle, y_circle_upper, 'r-', linewidth=2)
ax1.plot(x_circle, y_circle_lower, 'r-', linewidth=2)

ax1.set_xlim(-2, 2)
ax1.set_ylim(-0.5, 2.25)

#  Сетка с шагом 0.25 
ax1.set_xticks(np.arange(-2, 2, 0.25), minor=True)
ax1.set_yticks(np.arange(-0.5, 2.25, 0.25), minor=True)
# Подписи с шагом 0.5 
ax1.set_xticks(np.arange(-2, 2, 0.5))
ax1.set_yticks(np.arange(-0.5, 2.25, 0.5))
# Включаем сетку для обоих типов делений
ax1.grid(which='minor', alpha=0.5)  
ax1.grid(which='major', alpha=0.5) 

ax1.legend()
ax1.set_title('Графики функций')
ax1.set_xlabel('x')
ax1.set_ylabel('y')

ax2 = axes[1]
X, Y = np.meshgrid(np.linspace(-2, 2, 400), np.linspace(-0.5, 2.25, 400))

F1 = np.cos(X) - Y           
F2 = X**2 + (Y-1)**2 - 1     
residual = F1**2 + F2**2

# Рисуем линии уровня
colors = ['blue','green','lightyellow']
cmap_custom = LinearSegmentedColormap.from_list('yellow_red', colors, N=256)
im = ax2.pcolormesh(X, Y, residual, cmap=cmap_custom, shading='gouraud',norm=LogNorm(vmin=0.01, vmax=2))
plt.colorbar(im, ax=ax2, label='Невязки')

ax2.plot(x, y_cos, 'b--', alpha=0.4, linewidth=1)
ax2.plot(x_circle, y_circle_upper, 'r--', alpha=0.4, linewidth=1)
ax2.plot(x_circle, y_circle_lower, 'r--', alpha=0.4, linewidth=1)

ax2.set_title('Сумма квадратов невязок')
ax2.set_xlabel('x')
ax2.set_ylabel('y')

plt.tight_layout()
plt.savefig('2.png', dpi=150, bbox_inches='tight')
plt.show()