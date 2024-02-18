#include <iostream>
#include <vector>

using namespace std;

// Функция для вывода матрицы
void printMatrix(vector<vector<double>> &matrix)
{
    int n = matrix.size();
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            cout << matrix[i][j] << " ";
        }
        cout << endl;
    }
}

// Функция для создания единичной матрицы
vector<vector<double>> createIdentityMatrix(int n)
{
    vector<vector<double>> identity(n, vector<double>(n, 0));
    for (int i = 0; i < n; i++)
    {
        identity[i][i] = 1;
    }
    return identity;
}

// Функция для нахождения обратной матрицы методом Гаусса с выбором главного элемента по строке
vector<vector<double>> inverseMatrix(vector<vector<double>> &matrix)
{
    int n = matrix.size();

    // Создаем единичную матрицу
    vector<vector<double>> identity = createIdentityMatrix(n);

    for (int i = 0; i < n; i++)
    {
        // Поиск максимального элемента в столбце
        int maxRow = i;
        for (int j = i + 1; j < n; j++)
        {
            if (abs(matrix[j][i]) > abs(matrix[maxRow][i]))
            {
                maxRow = j;
            }
        }

        // Перестановка строк, если максимальный элемент не на диагонали
        if (maxRow != i)
        {
            swap(matrix[i], matrix[maxRow]);
            swap(identity[i], identity[maxRow]);
        }

        // Проверка того что главный элемент не равен 0, иначе определитель 0
        double pivot = matrix[i][i];
        if (pivot - 0 < 0.0001)
        {
            cout << "Определитель матрицы равен 0!!!\n";
            exit(1);
        }
        
        // Деление строки на главный элемент
        for (int j = i; j < n; j++)
        {
            matrix[i][j] /= pivot;
        }
        for (int j = 0; j < n; j++)
        {
            identity[i][j] /= pivot;
        }

        // Вычитание строк
        for (int k = 0; k < n; k++)
        {
            if (k != i)
            {
                double factor = matrix[k][i];
                for (int j = i; j < n; j++)
                {
                    matrix[k][j] -= factor * matrix[i][j];
                }
                for (int j = 0; j < n; j++)
                {
                    identity[k][j] -= factor * identity[i][j];
                }
            }
        }
    }

    return identity;
}

int main()
{
    int n;
    cout << "Введите размерность квадратной матрицы: ";
    cin >> n;

    cout << "Введите элементы матрицы:" << endl;
    vector<vector<double>> matrix(n, vector<double>(n));
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            cout << "[" << i << "]"
                 << "[" << j << "]:";
            cin >> matrix[i][j];
        }
    }

    vector<vector<double>> inverse = inverseMatrix(matrix);

    cout << "Обратная матрица:" << endl;
    printMatrix(inverse);

    return 0;
}
