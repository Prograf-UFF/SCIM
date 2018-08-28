#pragma once 

#include <vector>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <fstream>

namespace ds {

    /* Base class for dominance data.
     *
     * This class keeps a (N x n) dominance table E, where N is the number of subjects and n is the
     * number of stimuli.
     *
     * The element E(i,j) indicates the number of times subject i ranked stimulus j before other
     * stimuli minus the number of times the subject ranked it after other stimuli.
     */
    template<class MatrixType>
    class dominance_data {
    public:

        typedef MatrixType                   matrix_type;
        typedef const MatrixType&            matrix_const_reference;

        typedef typename matrix_type::Scalar real_type;
        typedef typename matrix_type::Index  size_type;

    public:

        virtual
        matrix_const_reference E() const = 0;
        
        inline
        size_type N() const {
            return E().rows();
        }

        inline
        size_type n() const {
            return E().cols();
        }
    };

    /* Rank-order data is a (N x n) matrix with the form:
     *
     *               x1 x2 x3 x4 x5 x6 x7 x8 x9
     *             +----------------------------+
     *   Subject 1 |  1  6  8  9  2  5  3  7  4 |
     *   Subject 2 |  6  9  8  5  3  1  7  2  4 |
     *   Subject 3 |  8  7  4  3  5  6  9  2  1 | Notice that the columns are ranked by the subjects.
     *             +----------------------------+
     *
     * where N is the number of subjects and n is the number of stimuli.
     */
    template<class MatrixType>
    class rank_order_data :
        public dominance_data<MatrixType> {
    private:

        typedef dominance_data<MatrixType> _super;

    public:

        typedef typename _super::matrix_type            matrix_type;
        typedef typename _super::matrix_const_reference matrix_const_reference;

    public:

        template<class Derived>
        inline
        rank_order_data(const Eigen::EigenBase<Derived> &data) :
            _super(),
            E_(data) {

            E_ *= -2;
            E_ += matrix_type::Constant(E_.rows() ,E_.cols(), E_.cols() + 1);
        }

    public:

        virtual
        matrix_const_reference E() const {
            return E_;
        }

    private:

        matrix_type E_;
    };

    /* Paired comparison data is defined by a collection of type std::vector<subject_opinions_type>,
     * where each element corresponds to the opinions of a subject.
     *
     * The opinions of a subject are defined by a collection of type std::vector<pair_type>, where
     * each element indicated whether stimulus pair_type::first was preferred over stimulus
     * pair_type::second.
     */
    template<class MatrixType>
    class paired_comparison_data :
        public dominance_data<MatrixType> {
    private:

        typedef dominance_data<MatrixType> _super;

    public:

        typedef typename _super::matrix_type            matrix_type;
        typedef typename _super::matrix_const_reference matrix_const_reference;

        typedef typename _super::size_type              size_type;

        typedef std::pair<size_type, size_type>         pair_type;
        typedef std::vector<pair_type>                  subject_opinions_type;
        typedef std::vector<subject_opinions_type>      subjects_opinions_type;

    public:

        inline
        paired_comparison_data(const subjects_opinions_type &data, size_type number_of_stimuli) :
            _super(),
            E_(data.size(), number_of_stimuli) {

            E_.setZero();

            for (size_t i1 = 0; i1 < data.size(); ++i1) {
                const subject_opinions_type &s = data[i1];

                for (size_t j = 0; j < s.size(); ++j) {
                    const pair_type &pair = s[j];

                    E_(i1, pair.first)++;
                    E_(i1, pair.second)--;
                }
            }
        }

    public:

        virtual
        matrix_const_reference E() const {
            return E_;
        }

    private:

        matrix_type E_;
    };

    /* Base class for incidence data.
     *
     * This class keeps a (n x m) incidence table F, where n is the number of groups and m the
     * number of stimuli.
     */
    template<class MatrixType>
    class incidence_data {
    public:

        typedef MatrixType                   matrix_type;
        typedef const MatrixType&            matrix_const_reference;

        typedef typename matrix_type::Scalar real_type;
        typedef typename matrix_type::Index  size_type;

    public:

        inline
        matrix_const_reference F() const {
            return F_;
        }
        
        inline
        size_type n() const {
            return F_.rows();
        }

        inline
        size_type m() const {
            return F_.cols();
        }

    protected:

        template<class Derived>
        inline
        incidence_data(const Eigen::EigenBase<Derived> &data) :
            F_(data) {
        }

    private:

        matrix_type F_;
    };

    /* A (n x m) contingency table has the form:
     *
     *               x1 x2 x3 x4 x5 x6 x7 x8 x9
     *             +----------------------------+
     *     Group 1 | 10 09 05 06 03 12 44 16 78 |
     *     Group 2 | 01 00 86 07 64 55 75 83 12 |
     *     Group 3 | 48 73 02 06 07 04 98 32 41 | -> 41 subjects in Group 3 choose x9.
     *             +----------------------------+
     *
     * where n is the number of groups and m the number of stimuli.
     */
    template<class MatrixType>
    class contingency_table :
        public incidence_data<MatrixType> {
    private:

        typedef incidence_data<MatrixType> _super;

    public:

        template<class Derived>
        inline
        contingency_table(const Eigen::EigenBase<Derived> &data) :
            _super(data) {
        }
    };

    /* Multiple-choice data is a (n x m) response-pattern matrix with the form:
     *
     *             +---------+---------+---------+
     *             |  Item 1 |  Item 2 |  Item 3 |
     *             +---------+---------+---------+
     *             | a  b  c | a  b  c | a  b  c |
     *             +---------+---------+---------+
     *   Subject 1 | 1  0  0 | 0  1  0 | 1  0  0 |
     *   Subject 2 | 0  1  0 | 1  0  0 | 0  0  1 |
     *   Subject 3 | 1  0  0 | 0  0  1 | 1  0  0 |
     *             +---------+---------+---------+
     *
     * where 1 means selection, n is the number of subjects and m is the number of possible
     * choices.
     */

    template<class MatrixType>
    class multiple_choice_data :
        public incidence_data<MatrixType> {
    private:

        typedef incidence_data<MatrixType> _super;

    public:

        template<class Derived>
        inline
        multiple_choice_data(const Eigen::EigenBase<Derived> &data) :
            _super(data) {
        }
    };

    /* Dual scaling procedure for dominance data.
     *
     * x_normed and y_normed are the normed weights for, respectively, the columns and the rows of
     * input data matrix with respect to each non-trivial solution.
     *
     * x_projected and y_projected are the weights for, respectively, the columns and the rows of
     * input data matrix with respect to each non-trivial solution.
     *
     * rho is the correlation ratio associated with each non-trivial solution.
     *
     * delta is the total variance explained by each solution.
     */
    template<class MatrixType, int Stimuli, int Dims, int XOptions, int MaxStimuli, int MaxDims, int Subjects, int YOptions, int MaxSubjects, int VectorOptions>
    void dual_scaling(
        const dominance_data<MatrixType> &data,
        Eigen::Matrix<typename MatrixType::Scalar, Stimuli, Dims, XOptions, MaxStimuli, MaxDims> &x_normed,
        Eigen::Matrix<typename MatrixType::Scalar, Subjects, Dims, YOptions, MaxSubjects, MaxDims> &y_normed,
        Eigen::Matrix<typename MatrixType::Scalar, Stimuli, Dims, XOptions, MaxStimuli, MaxDims> &x_projected,
        Eigen::Matrix<typename MatrixType::Scalar, Subjects, Dims, YOptions, MaxSubjects, MaxDims> &y_projected,
        Eigen::Matrix<typename MatrixType::Scalar, 1, Dims, VectorOptions, 1, MaxDims> &rho,
        Eigen::Matrix<typename MatrixType::Scalar, 1, Dims, VectorOptions, 1, MaxDims> &delta
    ) {
        typedef dominance_data<MatrixType>                                               dominance_data_type;

        typedef typename dominance_data_type::real_type                                  real_type;
        typedef typename dominance_data_type::size_type                                  size_type;
        
        typedef typename dominance_data_type::matrix_const_reference                     dominance_matrix_const_reference;

        typedef Eigen::Matrix<real_type, Stimuli, Dims, XOptions, MaxStimuli, MaxDims>   resulting_x_matrix_type;
        typedef Eigen::Matrix<real_type, 1, Dims, VectorOptions, 1, MaxDims>             resulting_vector_type;

        typedef Eigen::Matrix<real_type, 1, Stimuli>                                     row_vector_type;
        typedef typename row_vector_type::ConstantReturnType                             scalar_row_vector_type;

        typedef Eigen::Matrix<real_type, Stimuli, Stimuli>                               residual_matrix_type;

        typedef Eigen::SelfAdjointEigenSolver<residual_matrix_type>                      eigen_solver_type;
        typedef typename eigen_solver_type::MatrixType                                   eigenvectors_matrix_type;
        typedef typename eigen_solver_type::RealVectorType                               eigenvalues_vector_type;

        // Get arguments.
        const dominance_matrix_const_reference E = data.E();
        const size_type N = data.N();
        const size_type n = data.n();

        // Find eigenvalues and eigenvectors from generalized eigenequation of the residual matrix.
        residual_matrix_type C(E.transpose() * E);
        C /= (n * N * (n - 1) * (n - 1));

        const eigen_solver_type eigen_solver(C, Eigen::ComputeEigenvectors);
        const eigenvectors_matrix_type &V = eigen_solver.eigenvectors();
        const eigenvalues_vector_type &D = eigen_solver.eigenvalues(); // The eigenvalues are sorted in increasing order.

        // Get squared correlation ratio associated with the k-th non-trivial solution.
        rho.resize(n - 1);
        for (size_type i = 0; i < (n - 1); ++i) {
            rho[i] = sqrt(D[n - i - 1]);
        }
    
        // Get the weight of the columns.
        resulting_x_matrix_type x(n, n - 1);

        for (size_type i1 = 0; i1 < n; ++i1) {
            for (size_type i2 = 0; i2 < (n - 1); ++i2) {
                x(i1, i2) = -V(i1, n - i2 - 1);
            }
        }
    
        // Compute the constant multiplier for adjusting the unit of x.
        resulting_x_matrix_type x_sqr(n, n - 1);

        for (size_type i1 = 0; i1 < n; ++i1) {
            for (size_type i2 = 0; i2 < (n - 1); ++i2) {
                x_sqr(i1, i2) = x(i1, i2) * x(i1, i2);
            }
        }

        const real_type ft = N * n * (n - 1); // The total number of responses.

        const Eigen::DiagonalWrapper<const scalar_row_vector_type> Dc(row_vector_type::Constant(n, N * (n - 1))); // The marginal frequency of responses for cols.
        const resulting_x_matrix_type &T = Dc * x_sqr;

        resulting_vector_type cc(n - 1);

        for (size_type i2 = 0; i2 < (n - 1); ++i2) {
            real_type sum = 0;
            for (size_type i1 = 0; i1 < n; ++i1) {
                sum += T(i1, i2);
            }

            cc[i2] = sqrt(ft / sum);
        }

        // Compute the normed weights for columns.
        x_normed.resize(n, n - 1);

        for (size_type i1 = 0; i1 < n; ++i1) {
            for (size_type i2 = 0; i2 < (n - 1); ++i2) {
                x_normed(i1, i2) = x(i1, i2) * cc[i2];
            }
        }

        // Compute the normed weights for rows.
        const real_type fr = n * (n - 1); // The marginal frequency of responses for rows.

        y_normed = E * x_normed;

        for (size_type i2 = 0; i2 < (n - 1); ++i2) {
            const real_type t = 1 / (rho[i2] * fr);

            for (size_type i1 = 0; i1 < N; ++i1) {
                y_normed(i1, i2) *= t;
            }
        }

        // Compute the projected weights.
        x_projected.resize(n, n - 1);

        for (size_type i1 = 0; i1 < n; ++i1) {
            for (size_type i2 = 0; i2 < (n - 1); ++i2) {
                x_projected(i1, i2) = x_normed(i1, i2) * rho[i2];
            }
        }

        y_projected.resize(N, n - 1);

        for (size_type i1 = 0; i1 < N; ++i1) {
            for (size_type i2 = 0; i2 < (n - 1); ++i2) {
                y_projected(i1, i2) = y_normed(i1, i2) * rho[i2];
            }
        }

        // Compute the total variance explained by each solution.
        delta.resize(n - 1);

        real_type sum_rho_sqr = 0;
        for (size_type i = 0; i < (n - 1); ++i) {
            sum_rho_sqr += (delta[i] = D[n - i - 1]);
        }

        delta *= 100 / sum_rho_sqr;
    }

    /* Dual scaling procedure for incidence data.
     *
     * x_normed and y_normed are the normed weights for, respectively, the columns and the rows of
     * input data matrix with respect to each non-trivial solution.
     *
     * x_projected and y_projected are the weights for, respectively, the columns and the rows of
     * input data matrix with respect to each non-trivial solution.
     *
     * rho is the correlation ratio associated with each non-trivial solution.
     *
     * delta is the total variance explained by each solution.
     */
    template<class MatrixType, int Stimuli, int Dims, int XOptions, int MaxStimuli, int MaxDims, int Groups, int YOptions, int MaxGroups, int VectorOptions>
    void dual_scaling(
        const incidence_data<MatrixType> &data,
        Eigen::Matrix<typename MatrixType::Scalar, Stimuli, Dims, XOptions, MaxStimuli, MaxDims> &x_normed,
        Eigen::Matrix<typename MatrixType::Scalar, Groups, Dims, YOptions, MaxGroups, MaxDims> &y_normed,
        Eigen::Matrix<typename MatrixType::Scalar, Stimuli, Dims, XOptions, MaxStimuli, MaxDims> &x_projected,
        Eigen::Matrix<typename MatrixType::Scalar, Groups, Dims, YOptions, MaxGroups, MaxDims> &y_projected,
        Eigen::Matrix<typename MatrixType::Scalar, 1, Dims, VectorOptions, 1, MaxDims> &rho,
        Eigen::Matrix<typename MatrixType::Scalar, 1, Dims, VectorOptions, 1, MaxDims> &delta,
		Eigen::Matrix<typename MatrixType::Scalar, 1, Dims, VectorOptions, 1, MaxDims> &fc,
		Eigen::Matrix<typename MatrixType::Scalar, Groups, 1, YOptions, MaxGroups, 1> &fr,
		double LIMITE_INFERIOR_DELTA
    ) {
        typedef incidence_data<MatrixType>                                             incidence_data_type;

        typedef typename incidence_data_type::real_type                                real_type;
        typedef typename incidence_data_type::size_type                                size_type;
        
        typedef typename incidence_data_type::matrix_const_reference                   incidence_matrix_const_reference;

        typedef Eigen::Matrix<real_type, Stimuli, Dims, XOptions, MaxStimuli, MaxDims> resulting_x_matrix_type;
        typedef Eigen::Matrix<real_type, 1, Dims, VectorOptions, 1, MaxDims>           resulting_vector_type;

        typedef Eigen::Matrix<real_type, Groups, 1>                                    column_vector_type;
        typedef Eigen::Matrix<real_type, 1, Stimuli>                                   row_vector_type;

        typedef Eigen::Matrix<real_type, Stimuli, Stimuli>                             residual_matrix_type;

        // Get arguments.
        const incidence_matrix_const_reference F = data.F();
        size_type n = data.n();
        size_type m = data.m();

        // Compute marginal frequency of responses for rows and columns as (n x 1) and (m x 1) vectors, respectively, and the total number of responses.
        fr.resize(n);
        fc.resize(m);

        real_type ft = 0;

        for (size_type i1 = 0; i1 < n; ++i1) {
            fr[i1] = 0;
            for (size_type i2 = 0; i2 < m; ++i2) {
                fr[i1] += F(i1, i2);
                ft += F(i1, i2);
            }
        }

        for (size_type i2 = 0; i2 < m; ++i2) {
            fc[i2] = 0;
            for (size_type i1 = 0; i1 < n; ++i1) {
                fc[i2] += F(i1, i2);
            }

            if (fc[i2] == static_cast<real_type>(0)) {
                // The marginal frequence response of some columns is zero. This is a workaround.
                fc[i2] = std::numeric_limits<real_type>::epsilon();
                
                fr += column_vector_type::Constant(n, std::numeric_limits<real_type>::epsilon() / n);
                ft += std::numeric_limits<real_type>::epsilon();
            }
        }

        // Convert marginal frequencies of responses to diagonal matrices.
        const Eigen::DiagonalWrapper<column_vector_type> Dr(fr);
        const Eigen::DiagonalWrapper<row_vector_type> Dc(fc);

        residual_matrix_type M = F.transpose() * Dr.inverse() * F * Dc.inverse();
        Eigen::EigenSolver<residual_matrix_type> es(M);
        Eigen::MatrixXcd AD = es.eigenvalues();
        Eigen::MatrixXcd AV = es.eigenvectors();

        // sort autovector
		struct data {
			real_type number;
			size_t index;

			data(real_type n, size_t i){
				number = n;
				index = i;
			}
		};

		struct by_number {
			bool operator()(data const &left, data const &right) {
				return left.number > right.number;
			}
		};

		std::vector<data> sorted;
		
		for (size_type i = 0; i < AD.rows(); ++i)
			sorted.push_back(data(AD(i, 0).real(), i));

		std::sort(sorted.begin(), sorted.end(), by_number());

        //sort matrix of autovector
		Eigen::MatrixXcd VSort = Eigen::MatrixXcd(AV.rows(),AV.cols());
		int c = 0;
		for (data it : sorted) {
			Eigen::VectorXcd m_vec(AV.col(it.index));
			VSort.col(c++) = m_vec;
		}

        //remove irrelevant dimension
		int maxDimensao = 0;
        for (data it : sorted)
            if (it.number > LIMITE_INFERIOR_DELTA)
                 maxDimensao++;

        Eigen::MatrixXcd Vi = VSort.block(0, 1, AV.rows(), maxDimensao-1);
		residual_matrix_type V = Vi.real();
		Eigen::VectorXd D(maxDimensao-1);
		for (int i = 1; i < maxDimensao; ++i)
			D[i-1] = sorted[i].number;

		size_type dim = V.cols();

        // Get squared correlation ratio associated with the k-th non-trivial solution.
		rho.resize(dim);
		for (size_type i = 0; i < dim; ++i) {
			rho[i] = sqrt(D[i]);
        }

		 // Get the weight of the columns.
		resulting_x_matrix_type x(m, dim);

        for (size_type i1 = 0; i1 < m; ++i1) {
			for (size_type i2 = 0; i2 < dim; ++i2) {
                x(i1, i2) = V(i1, i2);
            }
        }

        // Compute the constant multiplier for adjusting the unit of x.
        resulting_x_matrix_type x_sqr(m, dim);

        for (size_type i1 = 0; i1 < m; ++i1) {
            for (size_type i2 = 0; i2 < dim; ++i2) {
                x_sqr(i1, i2) = x(i1, i2) * x(i1, i2);
            }
        }

        const resulting_x_matrix_type &T = Dc * x_sqr;


        resulting_vector_type cc(dim);

        for (size_type i2 = 0; i2 < dim; ++i2) {
            real_type sum = 0;
            for (size_type i1 = 0; i1 < m; ++i1) {
                sum += T(i1, i2);
            }

            cc[i2] = sqrt(ft / sum);
        }

        // Compute the normed weights for columns.
        x_normed.resize(m, dim);

        for (size_type i1 = 0; i1 < m; ++i1) {
            for (size_type i2 = 0; i2 < dim; ++i2) {
                x_normed(i1, i2) = x(i1, i2) * cc[i2];
            }
        }

        // Compute the normed weights for rows.
        y_normed = F * x_normed;

        for (size_type i2 = 0; i2 < dim; ++i2) {
            for (size_type i1 = 0; i1 < n; ++i1) {
                y_normed(i1, i2) /= (rho[i2] * fr[i1]);
            }
        }

        // Compute the projected weights.
        x_projected.resize(m, dim);

        for (size_type i1 = 0; i1 < m; ++i1) {
            for (size_type i2 = 0; i2 < dim; ++i2) {
                x_projected(i1, i2) = x_normed(i1, i2) * rho[i2];
            }
        }

        y_projected.resize(n, dim);

        for (size_type i1 = 0; i1 < n; ++i1) {
            for (size_type i2 = 0; i2 < dim; ++i2) {
                y_projected(i1, i2) = y_normed(i1, i2) * rho[i2];
            }
        }

        // Compute the total variance explained by each solution.
        delta.resize(dim);

        real_type sum_rho_sqr = 0;
        for (size_type i = 0; i < dim; ++i) {
            sum_rho_sqr += (delta[i] = D[i]);
        }

        delta *= 100 / sum_rho_sqr;
    }

}
