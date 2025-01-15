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
#include <geogram/mesh/mesh_geometry.h>

#include <assert.h>

using namespace GEO;

bool should_flip(vec3 p0, vec3 p1, vec3 p2, vec3 p3) {

    //p0-p1 is the edge 
    //p2 and p3 are the two opposite vertices 

    auto angle_between_three_vertices =
        [](vec3& S, vec3& M, vec3& Q) {

        vec3 p1 = S - M;
        vec3 p2 = Q - M;
        auto dot_pro = dot(p1, p2);

        return acos(dot_pro / (length(p1) * length(p2)));
        };


    // first check if the edge formed by v0-v1 is a delaunay edge
    // where v2 and v3 are the opposite vertices to the edge
    /*
        0
      / | \
     3  |  2
     \  |  /
        1
    */
    // if not delaunay, then we check if flipping it won't create a
    // foldover The case below would create a fold over
    /*
          0
        / | \
       /  1  \
      / /  \  \
      2       3
    */

    auto lambda = angle_between_three_vertices(p0, p2, p1);
    auto gamma = angle_between_three_vertices(p0, p3, p1);

    constexpr float PII = 3.14159265358979323f;

    if (lambda + gamma > PII + std::numeric_limits<float>::epsilon()) {
        // check if flipping won't create foldover

        auto alpha0 = angle_between_three_vertices(p3, p0, p1);
        auto beta0 = angle_between_three_vertices(p2, p0, p1);

        auto alpha1 = angle_between_three_vertices(p3, p1, p0);
        auto beta1 = angle_between_three_vertices(p2, p1, p0);

        if (alpha0 + beta0 < PII - std::numeric_limits<float>::epsilon() &&
            alpha1 + beta1 < PII - std::numeric_limits<float>::epsilon()) {
            return true;
        }
    }

    return false;
}


index_t find_thrid_vertex(Mesh& M, index_t f, index_t v1, index_t v2) {
    //given a face and two vertices of this face, find the third vertex of this face

    for (index_t c : M.facets.corners(f)) {
        if (M.facet_corners.vertex(c) != v1 && M.facet_corners.vertex(c) != v2) {
            return M.facet_corners.vertex(c);
        }
    }
    return NO_VERTEX;
}

bool flip_face_edges(Mesh& M) {
    //iterate over edges of this face and flip any if needed

    bool ret = false;

    for (index_t f : M.facets) {
        for (index_t c1 : M.facets.corners(f)) {
            index_t opposite_f = M.facet_corners.adjacent_facet(c1);
            if (opposite_f != NO_FACET) {
                index_t c2 = M.facets.next_corner_around_facet(f, c1);

                index_t c3 = M.facets.next_corner_around_facet(f, c2);

                //the edge vertices 
                index_t v1 = M.facet_corners.vertex(c1);
                index_t v2 = M.facet_corners.vertex(c2);

                //opposite vertices 
                index_t v3 = M.facet_corners.vertex(c3);
                index_t v4 = find_thrid_vertex(M, opposite_f, v1, v2);

                assert(v4 != NO_VERTEX);

                vec3 p1 = Geom::mesh_vertex(M, v1);
                vec3 p2 = Geom::mesh_vertex(M, v2);
                vec3 p3 = Geom::mesh_vertex(M, v3);
                vec3 p4 = Geom::mesh_vertex(M, v4);

                if (should_flip(p1, p2, p3, p4)) {

                    // TODO
                    // 
                    //M.facets.set_vertex(f, 0, v1);
                    //M.facets.set_vertex(f, 1, v4);
                    //M.facets.set_vertex(f, 2, v3);
                    //
                    //M.facets.set_vertex(opposite_f, 0, v4);
                    //M.facets.set_vertex(opposite_f, 1, v2);
                    //M.facets.set_vertex(opposite_f, 2, v3);

                    //also need to change "adjacent facets" so that "adjacent_facet" points to the right face next time
                    //  for(index_t c: M.facet_corners) {
                    //   M.facet_corners.set_adjacent_facet(c, FFF);
                    // }

                    // also need to change corner vertex 
                    //  for(index_t c: M.facet_corners) {
                    //   M.facet_corners.set_vertex(c, VVV);
                    // }
                    //

                    //ret = true;
                }


            }
        }
    }

    return ret;

}


int main(int argc, char** argv) {

    GEO::initialize(GEO::GEOGRAM_INSTALL_ALL);

    try {

        std::vector<std::string> filenames;

        if (!CmdLine::parse(argc, argv, filenames, "inmesh <outmesh>")) {
            return 1;
        }

        Mesh M;

        MeshIOFlags flags;
        flags.reset_element(MESH_CELLS);
        flags.set_attributes(MESH_ALL_ATTRIBUTES);
        if (!mesh_load(filenames[0], M, flags)) {
            return 1;
        }

        std::cout << "\n #### Input #Faces= " << M.facets.nb()
            << " #Vertices=" << M.vertices.nb() << "\n";

        auto start = std::chrono::high_resolution_clock::now();


        while (flip_face_edges(M)) {}

        auto stop = std::chrono::high_resolution_clock::now();

        std::cout << "\n $$$$ Delaunay Edge Flip took = "
            << std::chrono::duration<float, std::milli>(stop - start).count()
            << "(ms)\n";

        if (!mesh_save(M, filenames[1], flags)) {
            return 1;
        }

    }
    catch (const std::exception& e) {
        std::cerr << "Received an exception: " << e.what() << std::endl;
        return 1;
    }

    Logger::out("") << "Everything OK, Returning status 0" << std::endl;
    return 0;
}
