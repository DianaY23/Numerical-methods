#include <stdio.h>
#include <math.h>

// Структура для хранения результатов
typedef struct {
    double h;     // Шаг
    double error; // Ошибка
    double order; // Порядок сходимости
} Result;

// Функция для вычисления точного решения
double exact_solution(double lambda, double t) {
    return exp(lambda * t);
}

// Функция f(t, u) = λ * u
double func(double u, int lambda) {
    return lambda * u;
}

// Реализация метода Адамса второго порядка
void adams_second_order(double lambda, double h, int steps, Result *results, int index) {
    double time = 0.0;           // Время
    double u_prev = 1.0;         // Начальное условие
    double u_curr;               // Текущее значение
    double u_next;               // Следующее значение
    double exact_value;          // Точное решение
    double max_error = 0.0;      // Максимальная ошибка

    // Шаг методом Эйлера для получения первого значения
    u_curr = u_prev + h * lambda * u_prev;

    for (int step = 2; step <= steps; ++step) {
        time += h;

        // Явный метод Адамса второго порядка
        u_next = u_curr + h * lambda * (1.5 * u_curr - 0.5 * u_prev);

        // Вычисление ошибки
        exact_value = exact_solution(lambda, time);
        double error = fabs(u_next - exact_value);
        if (error > max_error) {
            max_error = error;
        }

        // Обновление значений для следующего шага
        u_prev = u_curr;
        u_curr = u_next;
    }

    // Сохранение результатов для текущего шага
    results[index].h = h;
    results[index].error = max_error;

    // Порядок сходимости будет вычислен позже
    results[index].order = -1.0; // Временное значение
}

// Реализация метода Рунге-Кутты
double runge_kutta(double h, int n, int lambda) {
    double u[n + 1], t[n + 1];
    u[0] = 1.0; // Начальное условие
    t[0] = 0.0;

    // Метод Рунге-Кутты
    for (int i = 1; i <= n; i++) {
        t[i] = t[i - 1] + h;
        u[i] = u[i - 1] * (1 + lambda * h / 2) / (1 - lambda * h / 2);
    }

    // Вычисление максимальной ошибки
    double max_error = fabs(u[0] - exp(lambda * t[0]));
    for (int i = 1; i < n; i++) {
        double error = fabs(u[i] - exp(lambda * t[i]));
        if (error > max_error) {
            max_error = error;
        }
    }
    
    return max_error;
}

// Функция для вычисления порядка сходимости для метода Адамса
void compute_orders(Result *results, int count) {
    for (int i = 1; i < count; ++i) {
        results[i].order = log(results[i - 1].error / results[i].error) / log(2.0);
    }
}

// Функция для вычисления порядка сходимости для метода Рунге-Кутты
void compute_runge_orders(double *h_values, int size, int lambda, Result *results) {
    for (int i = 0; i < size; i++) {
        double h = h_values[i];
        double error = runge_kutta(h, (int)(1 / h), lambda);
        results[i].h = h;
        results[i].error = error;

        if (i > 0) {
            results[i].order = log(results[i - 1].error / results[i].error) / log(2);
        } else {
            results[i].order = -1; // Для первого шага порядок не определён
        }
    }
}

// Функция для вывода результатов в таблицу
void print_results(const Result *results, int size) {
    printf("|    Step (h)    |    Error (ϕ)      |    Order (ψ)     |\n");
    printf("|----------------|-------------------|------------------|\n");
    for (int i = 0; i < size; ++i) {
        if (results[i].order >= 0) {
            printf("| %16.5f | %18.7f | %16.2f |\n",
                   results[i].h, results[i].error, results[i].order);
        } else {
            printf("| %16.5f | %18.7f |        ---        |\n",
                   results[i].h, results[i].error);
        }
    }
    printf("\n");
}

int main() {
    double lambda1 = -1.0;
    double lambda2 = -20.0;

    // Шаги для расчета (например 1/5, 1/10, 1/20, ...)
    double step_sizes[] = {0.2, 0.1, 0.05, 0.025, 0.0125, 0.00625};
    int num_sizes = sizeof(step_sizes) / sizeof(step_sizes[0]);

    Result adams_results1[6];
    Result adams_results2[6];
    Result runge_kutta_results1[6];
    Result runge_kutta_results2[6];

    // Метод Адамса второго порядка
    printf("Метод Адамса:\n");
    printf("Результаты для lambda = -1.0:\n");
    for (int i = 0; i < num_sizes; ++i) {
        int steps = (int)(2.0 / step_sizes[i]);
        adams_second_order(lambda1, step_sizes[i], steps, adams_results1, i);
    }
    compute_orders(adams_results1, num_sizes);
    print_results(adams_results1, num_sizes);

    printf("Результаты для lambda = -20.0:\n");
    for (int i = 0; i < num_sizes; ++i) {
        int steps = (int)(2.0 / step_sizes[i]);
        adams_second_order(lambda2, step_sizes[i], steps, adams_results2, i);
    }
    compute_orders(adams_results2, num_sizes);
    print_results(adams_results2, num_sizes);

    // Метод Рунге-Кутты
    printf("Метод Рунге-Кутты:\n");
    printf("Результаты для lambda = -1.0:\n");
    compute_runge_orders(step_sizes, num_sizes, -1, runge_kutta_results1);
    print_results(runge_kutta_results1, num_sizes);

    printf("Результаты для lambda = -20.0:\n");
    compute_runge_orders(step_sizes, num_sizes, -20, runge_kutta_results2);
    print_results(runge_kutta_results2, num_sizes);

    return 0;
}

