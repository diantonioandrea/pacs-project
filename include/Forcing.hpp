/**
 * @file Forcing.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-09
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef FORCING_PACS
#define FORCING_PACS

// Mesh.
#include <Mesh.hpp>

namespace pacs {

    // Function alias.
    using Function = double (*) (const double &, const double &);

    /**
     * @brief Source class for the Poisson problem.
     * 
     */
    class Source {
        private:
            
            // Function.
            Function function;

        public:

            // CONSTRUCTORS.

            Source(const Function &);

            // EVALUATION.

            double operator ()(const double &, const double &) const;
            Vector<double> operator ()(const Vector<double> &, const Vector<double> &) const;
    };

    // RHS.

    Vector<double> forcing(const Mesh &, const Source &);
}

#endif