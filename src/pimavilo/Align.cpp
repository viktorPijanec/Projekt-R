#include "Align.h"
#include <iostream>

namespace pimavilo
{

    int Align(
        const char *query, unsigned int query_len,
        const char *target, unsigned int target_len,
        AlignmentType type,
        int match,
        int mismatch,
        int gap,
        std::string *cigar,
        unsigned int *target_begin)
    {
        if (type == AlignmentType::Global)
        {
            // Creating a matrix
            Cell **matrix = new Cell *[query_len + 1];
            for (unsigned int i = 0; i <= query_len; ++i)
            {
                matrix[i] = new Cell[target_len + 1];
            }
            // Setting initial values
            matrix[0][0].score = 0;
            matrix[0][0].direction = Direction::Start;
            for (int i = 1; i <= query_len; ++i)
            {
                matrix[i][0].score = gap * i;
                matrix[i][0].direction = Direction::Left;
            }
            for (int i = 1; i <= target_len; ++i)
            {
                matrix[0][i].score = gap * i;
                matrix[0][i].direction = Direction::Up;
            }
            int match_mismatch;
            // Filling the matrix
            for (int i = 1; i <= query_len; ++i)
            {
                for (int j = 1; j <= target_len; ++j)
                {
                    if (query[i - 1] == target[j - 1])
                        match_mismatch = match;
                    else
                        match_mismatch = mismatch;

                    if (matrix[i - 1][j - 1].score + match_mismatch >= matrix[i - 1][j].score + gap &&
                        matrix[i - 1][j - 1].score + match_mismatch >= matrix[i][j - 1].score + gap)
                    {
                        matrix[i][j].score = matrix[i - 1][j - 1].score + match_mismatch;
                        matrix[i][j].direction = Direction::Diagonal;
                    }
                    else if (matrix[i - 1][j].score + gap >= matrix[i - 1][j - 1].score + match_mismatch &&
                             matrix[i - 1][j].score + match_mismatch >= matrix[i][j - 1].score + gap)
                    {
                        matrix[i][j].score = matrix[i - 1][j].score + gap;
                        matrix[i][j].direction = Direction::Up;
                    }
                    else
                    {
                        matrix[i][j].score = matrix[i][j - 1].score + gap;
                        matrix[i][j].direction = Direction::Left;
                    }
                }
            }

            int sol = matrix[query_len][target_len].score;
            // Free the allocated memory
            for (int i = 0; i <= query_len; ++i)
            {
                delete[] matrix[i];
            }
            delete[] matrix;
            return sol;
        }
        else if (type == AlignmentType::SemiGlobal)
        {
            // Creating a matrix
            Cell **matrix = new Cell *[query_len + 1];
            for (unsigned int i = 0; i <= query_len; ++i)
            {
                matrix[i] = new Cell[target_len + 1];
            }

            // Setting initial values
            matrix[0][0].score = 0;
            matrix[0][0].direction = Direction::Start;
            for (int i = 1; i <= query_len; ++i)
            {
                matrix[i][0].score = 0;
                matrix[i][0].direction = Direction::Start;
            }
            for (int i = 1; i <= target_len; ++i)
            {
                matrix[0][i].score = 0;
                matrix[0][i].direction = Direction::Start;
            }

            int match_mismatch;
            // Filling the matrix
            for (int i = 1; i <= query_len; ++i)
            {
                for (int j = 1; j <= target_len; ++j)
                {
                    if (query[i - 1] == target[j - 1])
                        match_mismatch = match;
                    else
                        match_mismatch = mismatch;

                    if (matrix[i - 1][j - 1].score + match_mismatch >= matrix[i - 1][j].score + gap &&
                        matrix[i - 1][j - 1].score + match_mismatch >= matrix[i][j - 1].score + gap)
                    {
                        matrix[i][j].score = matrix[i - 1][j - 1].score + match_mismatch;
                        matrix[i][j].direction = Direction::Diagonal;
                    }
                    else if (matrix[i - 1][j].score + gap >= matrix[i - 1][j - 1].score + match_mismatch &&
                             matrix[i - 1][j].score + match_mismatch >= matrix[i][j - 1].score + gap)
                    {
                        matrix[i][j].score = matrix[i - 1][j].score + gap;
                        matrix[i][j].direction = Direction::Up;
                    }
                    else
                    {
                        matrix[i][j].score = matrix[i][j - 1].score + gap;
                        matrix[i][j].direction = Direction::Left;
                    }
                }
            }

            // Computing the solution
            int sol = matrix[query_len][target_len].score;
            for (int i = 0; i <= target_len; ++i)
            {
                if (matrix[query_len][i].score >= sol)
                {
                    sol = matrix[query_len][i].score;
                }
            }
            for (int i = 0; i <= query_len; ++i)
            {
                if (matrix[i][target_len].score >= sol)
                {
                    sol = matrix[i][target_len].score;
                }
            }
            // Free the allocated memory
            for (int i = 0; i <= query_len; ++i)
            {
                delete[] matrix[i];
            }
            delete[] matrix;
            return sol;
        }
        else if (type == AlignmentType::Local)
        {
            // Creating a matrix
            Cell **matrix = new Cell *[query_len + 1];
            for (unsigned int i = 0; i <= query_len; ++i)
            {
                matrix[i] = new Cell[target_len + 1];
            }

            // Setting initial values
            matrix[0][0].score = 0;
            matrix[0][0].direction = Direction::Start;
            for (int i = 1; i <= query_len; ++i)
            {
                matrix[i][0].score = 0;
                matrix[i][0].direction = Direction::Start;
            }
            for (int i = 1; i <= target_len; ++i)
            {
                matrix[0][i].score = 0;
                matrix[0][i].direction = Direction::Start;
            }

            int match_mismatch;
            // Filling the matrix
            for (int i = 1; i <= query_len; ++i)
            {
                for (int j = 1; j <= target_len; ++j)
                {
                    if (query[i - 1] == target[j - 1])
                        match_mismatch = match;
                    else
                        match_mismatch = mismatch;

                    if (0 >= matrix[i - 1][j - 1].score + match_mismatch &&
                        0 >= matrix[i - 1][j].score + gap &&
                        0 >= matrix[i][j - 1].score + gap)
                    {
                        matrix[i][j].score = 0;
                        matrix[i][j].direction = Direction::Start;
                    }
                    else if (matrix[i - 1][j - 1].score + match_mismatch >= matrix[i - 1][j].score + gap &&
                             matrix[i - 1][j - 1].score + match_mismatch >= matrix[i][j - 1].score + gap)
                    {
                        matrix[i][j].score = matrix[i - 1][j - 1].score + match_mismatch;
                        matrix[i][j].direction = Direction::Diagonal;
                    }
                    else if (matrix[i - 1][j].score + gap >= matrix[i - 1][j - 1].score + match_mismatch &&
                             matrix[i - 1][j].score + match_mismatch >= matrix[i][j - 1].score + gap)
                    {
                        matrix[i][j].score = matrix[i - 1][j].score + gap;
                        matrix[i][j].direction = Direction::Up;
                    }
                    else
                    {
                        matrix[i][j].score = matrix[i][j - 1].score + gap;
                        matrix[i][j].direction = Direction::Left;
                    }
                }
            }

            // Computing the solution
            int sol = matrix[query_len][target_len].score;
            for (int i = 0; i <= query_len; ++i)
            {
                for (int j = 0; j <= target_len; ++j)
                {
                    if (matrix[i][j].score >= sol)
                    {
                        sol = matrix[i][j].score;
                    }
                }
            }
            // Free the allocated memory
            for (int i = 0; i <= query_len; ++i)
            {
                delete[] matrix[i];
            }
            delete[] matrix;
            return sol;
        }
        return 0;
    }

} // namespace pimavilo
