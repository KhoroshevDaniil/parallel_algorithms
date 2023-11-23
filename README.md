# Параллельные алгоритмы (лабы 2023)

**Цель лабораторных** - реализовать параллельное вычисление разностной схемы (в моем случае - неявной модифицированной) с помощью технологий OpenMP и MPI, а также векторизации последовательного алгоритма

# Описание задачи (разностной схемы)
Исходный отчёт можно посмотреть вот тут: [отчёт](/kursach.pdf)

На вход слабопоглощающего оптического элемента, представляющего
собой цилиндрический фрагмент оптоволокна кругового сечения длиной $l$ и
радиусом $R$, подаётся лазерное излучение с распределением интенсивности
$I_0(r)$. Ослабление (затухание) монохроматического лазерного пучка при его
распределении в поглощающей среде описывается законом Бугера
$$I(z, r) = I_0(r) exp(−\beta z),$$
где $I(z, r)$ – интенсивность излучения, прошедшего через слой вещества
толщиной $z$,

$\beta$ – коэффициент поглощения энергии излучения.

В результате поглощения части энергии излучения оптический элемент
нагревается.

Оптический элемент выполнен из однородного материала,
характеризуемого коэффициентом поглощения $\beta$, теплопроводности $k̃$,
объёмной теплоёмкости $c$. Пучок света обладает круговой симметрией и
падает нормально на входную грань оптоволокна, причём оси пучка и
оптического элемента совпадают.

Между поверхностями $r = R$, $z = 0$ и $z = l$ оптического элемента и
окружающей средой имеет место теплообмен, описываемый законом Ньютона
с коэффициентом теплообмена $\alpha$.
В момент включения лазера $t = 0$ температура оптического элемента
предполагается одинаковой во всех точках и равной температуре окружающей
среды $u_0$.

Разработать программу расчёта среднего по радиусу значения
температуры оптического элемента на временном промежутке $0 < t ⩽ T$.
# Компиляция и запуск программ
Для последовательного алгоритма
```
g++ -o cpu main.cpp scheme.cpp utils.cpp
./cpu
```
Оптимизация (в том числе - векторизация) последовательного алгоритма (O - буква, 3 - цифра)
```
g++ -O3 -o cpu_v main.cpp scheme.cpp utils.cpp
./cpu_v
```
Для OpenMP
```
g++ -fopenmp -o omp omp.cpp utils.cpp
./omp
```
Для MPI
```
mpic++ -o mpi mpi.cpp utils.cpp
mpirun -n 2 ./mpi
```
# Выводы
Время работы каждой программы было усреднено по 12 запускам для 4х пар параметров $K$ и $I$: $[200, 200], [800, 800], [1500, 1500], [3000, 3000]$

Под конец в лидеры по ускорению стала выбиваться OpenMP.

![графики времени и ускорения работы параллельных программ](/summary.jpeg)
