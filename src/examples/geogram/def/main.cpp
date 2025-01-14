/*
 *  Copyright (c) 2000-2022 Inria
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice,
 *  this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *  this list of conditions and the following disclaimer in the documentation
 *  and/or other materials provided with the distribution.
 *  * Neither the name of the ALICE Project-Team nor the names of its
 *  contributors may be used to endorse or promote products derived from this
 *  software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 *  Contact: Bruno Levy
 *
 *     https://www.inria.fr/fr/bruno-levy
 *
 *     Inria,
 *     Domaine de Voluceau,
 *     78150 Le Chesnay - Rocquencourt
 *     FRANCE
 *
 */

#include <geogram/basic/common.h>
#include <geogram/basic/logger.h>
#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>
#include <geogram/basic/stopwatch.h>
#include <geogram/basic/file_system.h>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_subdivision.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/NL/nl.h>

namespace {
    using namespace GEO;


    void mesh_smooth(Mesh* M, NLenum solver = NL_SOLVER_DEFAULT) {

        // Chain corners around vertices
        vector<index_t> v2c(M->vertices.nb(), index_t(-1));
        vector<index_t> next_c_around_v(M->facet_corners.nb(), index_t(-1));
        vector<index_t> c2f(M->facet_corners.nb(), index_t(-1));
        for(index_t f=0; f<M->facets.nb(); ++f) {
            for(index_t c=M->facets.corners_begin(f);
                c<M->facets.corners_end(f); ++c
               ) {
                index_t v = M->facet_corners.vertex(c);
                next_c_around_v[c] = v2c[v];
                v2c[v] = c;
                c2f[c] = f;
            }
        }

        nlNewContext();

        if(solver == NL_SUPERLU_EXT || solver == NL_PERM_SUPERLU_EXT) {
            if(nlInitExtension("SUPERLU")) {
                nlSolverParameteri(NL_SOLVER, NLint(solver));
            } else {
                Logger::err("NL") << "Could not init SUPERLU extension";
            }
        } else if(solver == NL_CHOLMOD_EXT) {
            if(nlInitExtension("CHOLMOD")) {
                nlSolverParameteri(NL_SOLVER, NLint(solver));
            } else {
                Logger::err("NL") << "Could not init CHOLMOD extension";
            }
        }

        nlSolverParameteri(NL_LEAST_SQUARES, NL_TRUE);
        nlSolverParameteri(NL_NB_VARIABLES, NLint(M->vertices.nb()));
        nlSolverParameteri(NL_NB_SYSTEMS, NLint(M->vertices.dimension()));
        nlEnable(NL_NORMALIZE_ROWS);
        nlEnable(NL_VARIABLES_BUFFER);

        Attribute<bool> v_is_locked(M->vertices.attributes(), "selection");

        nlBegin(NL_SYSTEM);

        for(index_t coord=0; coord<M->vertices.dimension(); ++coord) {
            // Bind directly the variables buffer to the coordinates in
            // the mesh, to avoid copying data.
            nlBindBuffer(
                NL_VARIABLES_BUFFER, NLuint(coord),
                M->vertices.point_ptr(0) + coord,
                NLuint(sizeof(double)*M->vertices.dimension())
            );
        }

        for(index_t v=0; v<M->vertices.nb(); ++v) {
            if(v_is_locked[v]) {
                nlLockVariable(v);
            }
        }

        nlBegin(NL_MATRIX);
        for(index_t v=0; v<M->vertices.nb(); ++v) {
            nlBegin(NL_ROW);
            index_t count = 0;
            for(
                index_t c = v2c[v];
                c != index_t(-1); c = next_c_around_v[c]
            ) {
                index_t f = c2f[c];
                index_t c2 = M->facets.next_corner_around_facet(f,c);
                index_t w = M->facet_corners.vertex(c2);
                nlCoefficient(w, 1.0);
                ++count;
            }
            nlCoefficient(v, -double(count));
            nlEnd(NL_ROW);
        }
        nlEnd(NL_MATRIX);
        nlEnd(NL_SYSTEM);

        Logger::div("Solve");
        nlSolve();

        nlDeleteContext(nlGetCurrent());
    }
}

/*
 * Turns out M->edges is not actually populated
 * during mesh construction!
 */
void delete_an_edges_test(Mesh* M) {
    vector<index_t> delete_e(M->edges.nb(),0);
    for(index_t e=0; e<M->edges.nb()/2; ++e) {
        delete_e[e] = 1;
    }
    M->edges.delete_elements(delete_e);
}

void delete_faces_test(Mesh& M) {
    vector<index_t> delete_f(M.facets.nb(),0);
    for(index_t f=0; f<M.facets.nb()/4; ++f) {
        delete_f[f] = 1;
    }
    M.facets.delete_elements(delete_f);
}

// just a work in progress
void angle_test() {
    // find the angle between S, M, Q vertices (i.e., angle at M)
    auto angle_between_three_vertices =
        [](vec3& S,
           vec3& M,
           vec3& Q) {

            vec3 p1      = S - M;
            vec3 p2      = Q - M;
            double dot_pro = dot(p1, p2);// p1.dot(p2);

            return acos(dot_pro / (length(p1) * length(p2)) );

            //        if constexpr (std::is_same_v<T, float>) {
            //            return acosf(dot_pro / (p1.norm() * p2.norm()));
            //        } else {
            //            return acos(dot_pro / (p1.norm() * p2.norm()));
            //        }
    };
    vec3 p0, p1, p2;
    auto angle = angle_between_three_vertices(p0, p1, p2);
}
void explore_mesh_connectivity(Mesh& M) {

    //for(index_t f=0; f<M.facets.nb(); ++f) {
    for (index_t f : M.facets) {
        std::cout << "f " << f << "\n";
        std::cout << "\tvids: ";
        // Look at a face's vertex indicies
        for (index_t lv = 0; lv < M.facets.nb_vertices(f); ++lv) {
            index_t v = M.facets.vertex(f, lv);
            std::cout << v << " ";
        }
        std::cout << "\n\n";

        // Look at a face's corners
        for (index_t c = M.facets.corners_begin(f); c != M.facets.corners_end(f); ++c) {
            std::cout << "\tcorner " << c
                << ": vertex " << M.facet_corners.vertex(c)
                << "\n";
        }
        std::cout << "\n";
        std::cout << "\tGiven a face's local edge, get an adjacent face.\n";
        std::cout << "\tThen, ask that face for a face adjacent to it.\n";
        // Given a local edge index, find a face adjacent to the edge
        for (index_t e = 0; e < 3; ++e) {
            auto fa = M.facets.adjacent(f, e);
            std::cout << "\te" << e << " fa: " << fa;
            std::cout << " fa " << ": " << M.facets.adjacent(fa, e) << std::endl;
        }

        for(index_t c = M.facets.corners_begin(f);
            c < M.facets.corners_end(f); ++c) {

            }
        std::cout << "----------" << std::endl;

        // std::cout << "f: " << f << " verts: " << M.facets.nb_vertices(f);
        // std::cout << " corners: " << M.facets.nb_corners(f) << std::endl;
    }

    // Output some mesh stats at the end here
    std::cout << "vertices.nb(): " << M.vertices.nb() << std::endl;
    std::cout << "facets.nb(): " << M.facets.nb() << std::endl;
    std::cout << "edges.nb(): " << M.edges.nb() << std::endl;
    std::cout << "facet_corners.nb(): " << M.facet_corners.nb() << std::endl;
}
int main(int argc, char** argv) {
    using namespace GEO;

    GEO::initialize(GEO::GEOGRAM_INSTALL_ALL);

    try {
        Stopwatch Wtot("Total time");

        std::vector<std::string> filenames;

        CmdLine::import_arg_group("standard");
        CmdLine::import_arg_group("algo");
        CmdLine::declare_arg(
            "solver", "NL_SOLVER_DEFAULT", "solver"
        );
        CmdLine::declare_arg(
            "nb_subdivide", 0, "number of times mesh is subdivided"
        );

        if(
            !CmdLine::parse(
                argc, argv, filenames, "inmesh <outmesh>"
            )
        ) {
            return 1;
        }


        if(filenames.size() != 2) {
            Logger::out("Smooth") << "Generating output to out.geogram"
                                  << std::endl;
            filenames.push_back("out.geogram");
        }

        Logger::div("Data I/O");

        Mesh M;

        MeshIOFlags flags;
        flags.reset_element(MESH_CELLS);
        flags.set_attributes(MESH_ALL_ATTRIBUTES);
        if(!mesh_load(filenames[0], M, flags)) {
            return 1;
        }

        explore_mesh_connectivity(M);

        std::cout << "delete some edges\n";
        delete_an_edges_test(&M);

        std::cout << "delete some faces\n";
        delete_faces_test(M);

        // {
        //     Attribute<bool> is_locked(M.vertices.attributes(), "selection");
        //     for(index_t v=0; v<M.vertices.nb(); ++v) {
        //         is_locked[v] = true;
        //     }
        //     for(index_t i=0; i<CmdLine::get_arg_uint("nb_subdivide"); ++i) {
        //         mesh_split_triangles(M);
        //     }
        // }

        // NLenum solver = NL_SOLVER_DEFAULT;
        // std::string solver_string = CmdLine::get_arg("solver");
        // if(solver_string == "NL_CG") {
        //     solver = NL_CG;
        // } else if(solver_string == "NL_SUPERLU_EXT") {
        //     solver = NL_SUPERLU_EXT;
        // } else if(solver_string == "NL_PERM_SUPERLU_EXT") {
        //     solver = NL_PERM_SUPERLU_EXT;
        // } else if(solver_string == "NL_SYMMETRIC_SUPERLU_EXT") {
        //     solver = NL_SYMMETRIC_SUPERLU_EXT;
        // } else if(solver_string == "NL_CHOLMOD_EXT") {
        //     solver = NL_CHOLMOD_EXT;
        // }

        // mesh_smooth(&M, solver);

        if(!mesh_save(M, filenames[1], flags)) {
            return 1;
        }

    }
    catch(const std::exception& e) {
        std::cerr << "Received an exception: " << e.what() << std::endl;
        return 1;
    }

    Logger::out("") << "Everything OK, Returning status 0" << std::endl;
    return 0;
}
