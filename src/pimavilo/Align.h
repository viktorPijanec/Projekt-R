#ifndef PIMAVILO_ALIGN_H
#define PIMAVILO_ALIGN_H

#include <string>

namespace pimavilo
{

    enum class AlignmentType
    {
        Global,
        Local,
        SemiGlobal
    };

    enum class Direction
    {
        Up,
        Left,
        Diagonal,
        Start
    };

    class Cell
    {
    public:
        int score;
        Direction direction;
    };

    int Align(
        const char *query, unsigned int query_len,
        const char *target, unsigned int target_len,
        AlignmentType type,
        int match,
        int mismatch,
        int gap,
        std::string *cigar = nullptr,
        unsigned int *target_begin = nullptr);

} // namespace pimavilo

#endif // PIMAVILO_ALIGN_H
