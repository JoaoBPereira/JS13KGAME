// Represents a 4x4 matrix stored in row-major order that uses Float32Arrays
// when available. Matrix operations can either be done using convenient
// methods that return a new matrix for the result or optimized methods
// that store the result in an existing matrix to avoid generating garbage.

var hasFloat32Array = (typeof Float32Array !== 'undefined');

// ### new GL.Matrix([elements])
//
// This constructor takes 16 arguments in row-major order, which can be passed
// individually, as a list, or even as four lists, one for each row. If the
// arguments are omitted then the identity matrix is constructed instead.
function Matrix() {
    var m = Array.prototype.concat.apply([], arguments);
    if (!m.length) {
        m = [
            1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 1
        ];
    }
    this.m = hasFloat32Array ? new Float32Array(m) : m;
}

Matrix.prototype = {
    // ### .inverse()
    //
    // Returns the matrix that when multiplied with this matrix results in the
    // identity matrix.
    inverse: function () {
        return Matrix.inverse(this, new Matrix());
    },

    // ### .transpose()
    //
    // Returns this matrix, exchanging columns for rows.
    transpose: function () {
        return Matrix.transpose(this, new Matrix());
    },

    // ### .multiply(matrix)
    //
    // Returns the concatenation of the transforms for this matrix and `matrix`.
    // This emulates the OpenGL function `glMultMatrix()`.
    multiply: function (matrix) {
        return Matrix.multiply(this, matrix, new Matrix());
    },

    // ### .transformPoint(point)
    //
    // Transforms the vector as a point with a w coordinate of 1. This
    // means translations will have an effect, for example.
    transformPoint: function (v) {
        var m = this.m;
        return new Vector(
            m[0] * v.x + m[1] * v.y + m[2] * v.z + m[3],
            m[4] * v.x + m[5] * v.y + m[6] * v.z + m[7],
            m[8] * v.x + m[9] * v.y + m[10] * v.z + m[11]
        ).divide(m[12] * v.x + m[13] * v.y + m[14] * v.z + m[15]);
    },

    // ### .transformPoint(vector)
    //
    // Transforms the vector as a vector with a w coordinate of 0. This
    // means translations will have no effect, for example.
    transformVector: function (v) {
        var m = this.m;
        return new Vector(
            m[0] * v.x + m[1] * v.y + m[2] * v.z,
            m[4] * v.x + m[5] * v.y + m[6] * v.z,
            m[8] * v.x + m[9] * v.y + m[10] * v.z
        );
    }
};

// ### GL.Matrix.inverse(matrix[, result])
//
// Returns the matrix that when multiplied with `matrix` results in the
// identity matrix. You can optionally pass an existing matrix in `result`
// to avoid allocating a new matrix. This implementation is from the Mesa
// OpenGL function `__gluInvertMatrixd()` found in `project.c`.
Matrix.inverse = function (matrix, result) {
    result = result || new Matrix();
    var m = matrix.m, r = result.m;

    r[0] = m[5] * m[10] * m[15] - m[5] * m[14] * m[11] - m[6] * m[9] * m[15] + m[6] * m[13] * m[11] + m[7] * m[9] * m[14] - m[7] * m[13] * m[10];
    r[1] = -m[1] * m[10] * m[15] + m[1] * m[14] * m[11] + m[2] * m[9] * m[15] - m[2] * m[13] * m[11] - m[3] * m[9] * m[14] + m[3] * m[13] * m[10];
    r[2] = m[1] * m[6] * m[15] - m[1] * m[14] * m[7] - m[2] * m[5] * m[15] + m[2] * m[13] * m[7] + m[3] * m[5] * m[14] - m[3] * m[13] * m[6];
    r[3] = -m[1] * m[6] * m[11] + m[1] * m[10] * m[7] + m[2] * m[5] * m[11] - m[2] * m[9] * m[7] - m[3] * m[5] * m[10] + m[3] * m[9] * m[6];

    r[4] = -m[4] * m[10] * m[15] + m[4] * m[14] * m[11] + m[6] * m[8] * m[15] - m[6] * m[12] * m[11] - m[7] * m[8] * m[14] + m[7] * m[12] * m[10];
    r[5] = m[0] * m[10] * m[15] - m[0] * m[14] * m[11] - m[2] * m[8] * m[15] + m[2] * m[12] * m[11] + m[3] * m[8] * m[14] - m[3] * m[12] * m[10];
    r[6] = -m[0] * m[6] * m[15] + m[0] * m[14] * m[7] + m[2] * m[4] * m[15] - m[2] * m[12] * m[7] - m[3] * m[4] * m[14] + m[3] * m[12] * m[6];
    r[7] = m[0] * m[6] * m[11] - m[0] * m[10] * m[7] - m[2] * m[4] * m[11] + m[2] * m[8] * m[7] + m[3] * m[4] * m[10] - m[3] * m[8] * m[6];

    r[8] = m[4] * m[9] * m[15] - m[4] * m[13] * m[11] - m[5] * m[8] * m[15] + m[5] * m[12] * m[11] + m[7] * m[8] * m[13] - m[7] * m[12] * m[9];
    r[9] = -m[0] * m[9] * m[15] + m[0] * m[13] * m[11] + m[1] * m[8] * m[15] - m[1] * m[12] * m[11] - m[3] * m[8] * m[13] + m[3] * m[12] * m[9];
    r[10] = m[0] * m[5] * m[15] - m[0] * m[13] * m[7] - m[1] * m[4] * m[15] + m[1] * m[12] * m[7] + m[3] * m[4] * m[13] - m[3] * m[12] * m[5];
    r[11] = -m[0] * m[5] * m[11] + m[0] * m[9] * m[7] + m[1] * m[4] * m[11] - m[1] * m[8] * m[7] - m[3] * m[4] * m[9] + m[3] * m[8] * m[5];

    r[12] = -m[4] * m[9] * m[14] + m[4] * m[13] * m[10] + m[5] * m[8] * m[14] - m[5] * m[12] * m[10] - m[6] * m[8] * m[13] + m[6] * m[12] * m[9];
    r[13] = m[0] * m[9] * m[14] - m[0] * m[13] * m[10] - m[1] * m[8] * m[14] + m[1] * m[12] * m[10] + m[2] * m[8] * m[13] - m[2] * m[12] * m[9];
    r[14] = -m[0] * m[5] * m[14] + m[0] * m[13] * m[6] + m[1] * m[4] * m[14] - m[1] * m[12] * m[6] - m[2] * m[4] * m[13] + m[2] * m[12] * m[5];
    r[15] = m[0] * m[5] * m[10] - m[0] * m[9] * m[6] - m[1] * m[4] * m[10] + m[1] * m[8] * m[6] + m[2] * m[4] * m[9] - m[2] * m[8] * m[5];

    var det = m[0] * r[0] + m[1] * r[4] + m[2] * r[8] + m[3] * r[12];
    for (var i = 0; i < 16; i++) r[i] /= det;
    return result;
};

// ### GL.Matrix.transpose(matrix[, result])
//
// Returns `matrix`, exchanging columns for rows. You can optionally pass an
// existing matrix in `result` to avoid allocating a new matrix.
Matrix.transpose = function (matrix, result) {
    result = result || new Matrix();
    var m = matrix.m, r = result.m;
    r[0] = m[0]; r[1] = m[4]; r[2] = m[8]; r[3] = m[12];
    r[4] = m[1]; r[5] = m[5]; r[6] = m[9]; r[7] = m[13];
    r[8] = m[2]; r[9] = m[6]; r[10] = m[10]; r[11] = m[14];
    r[12] = m[3]; r[13] = m[7]; r[14] = m[11]; r[15] = m[15];
    return result;
};

// ### GL.Matrix.multiply(left, right[, result])
//
// Returns the concatenation of the transforms for `left` and `right`. You can
// optionally pass an existing matrix in `result` to avoid allocating a new
// matrix. This emulates the OpenGL function `glMultMatrix()`.
Matrix.multiply = function (left, right, result) {
    result = result || new Matrix();
    var a = left.m, b = right.m, r = result.m;

    r[0] = a[0] * b[0] + a[1] * b[4] + a[2] * b[8] + a[3] * b[12];
    r[1] = a[0] * b[1] + a[1] * b[5] + a[2] * b[9] + a[3] * b[13];
    r[2] = a[0] * b[2] + a[1] * b[6] + a[2] * b[10] + a[3] * b[14];
    r[3] = a[0] * b[3] + a[1] * b[7] + a[2] * b[11] + a[3] * b[15];

    r[4] = a[4] * b[0] + a[5] * b[4] + a[6] * b[8] + a[7] * b[12];
    r[5] = a[4] * b[1] + a[5] * b[5] + a[6] * b[9] + a[7] * b[13];
    r[6] = a[4] * b[2] + a[5] * b[6] + a[6] * b[10] + a[7] * b[14];
    r[7] = a[4] * b[3] + a[5] * b[7] + a[6] * b[11] + a[7] * b[15];

    r[8] = a[8] * b[0] + a[9] * b[4] + a[10] * b[8] + a[11] * b[12];
    r[9] = a[8] * b[1] + a[9] * b[5] + a[10] * b[9] + a[11] * b[13];
    r[10] = a[8] * b[2] + a[9] * b[6] + a[10] * b[10] + a[11] * b[14];
    r[11] = a[8] * b[3] + a[9] * b[7] + a[10] * b[11] + a[11] * b[15];

    r[12] = a[12] * b[0] + a[13] * b[4] + a[14] * b[8] + a[15] * b[12];
    r[13] = a[12] * b[1] + a[13] * b[5] + a[14] * b[9] + a[15] * b[13];
    r[14] = a[12] * b[2] + a[13] * b[6] + a[14] * b[10] + a[15] * b[14];
    r[15] = a[12] * b[3] + a[13] * b[7] + a[14] * b[11] + a[15] * b[15];

    return result;
};

// ### GL.Matrix.identity([result])
//
// Returns an identity matrix. You can optionally pass an existing matrix in
// `result` to avoid allocating a new matrix. This emulates the OpenGL function
// `glLoadIdentity()`.
Matrix.identity = function (result) {
    result = result || new Matrix();
    var m = result.m;
    m[0] = m[5] = m[10] = m[15] = 1;
    m[1] = m[2] = m[3] = m[4] = m[6] = m[7] = m[8] = m[9] = m[11] = m[12] = m[13] = m[14] = 0;
    return result;
};

// ### GL.Matrix.perspective(fov, aspect, near, far[, result])
//
// Returns a perspective transform matrix, which makes far away objects appear
// smaller than nearby objects. The `aspect` argument should be the width
// divided by the height of your viewport and `fov` is the top-to-bottom angle
// of the field of view in degrees. You can optionally pass an existing matrix
// in `result` to avoid allocating a new matrix. This emulates the OpenGL
// function `gluPerspective()`.
Matrix.perspective = function (fov, aspect, near, far, result) {
    var y = Math.tan(fov * Math.PI / 360) * near;
    var x = y * aspect;
    return Matrix.frustum(-x, x, -y, y, near, far, result);
};

// ### GL.Matrix.frustum(left, right, bottom, top, near, far[, result])
//
// Sets up a viewing frustum, which is shaped like a truncated pyramid with the
// camera where the point of the pyramid would be. You can optionally pass an
// existing matrix in `result` to avoid allocating a new matrix. This emulates
// the OpenGL function `glFrustum()`.
Matrix.frustum = function (l, r, b, t, n, f, result) {
    result = result || new Matrix();
    var m = result.m;

    m[0] = 2 * n / (r - l);
    m[1] = 0;
    m[2] = (r + l) / (r - l);
    m[3] = 0;

    m[4] = 0;
    m[5] = 2 * n / (t - b);
    m[6] = (t + b) / (t - b);
    m[7] = 0;

    m[8] = 0;
    m[9] = 0;
    m[10] = -(f + n) / (f - n);
    m[11] = -2 * f * n / (f - n);

    m[12] = 0;
    m[13] = 0;
    m[14] = -1;
    m[15] = 0;

    return result;
};


// ### GL.Matrix.scale(x, y, z[, result])
//
// This emulates the OpenGL function `glScale()`. You can optionally pass an
// existing matrix in `result` to avoid allocating a new matrix.
Matrix.scale = function (x, y, z, result) {
    result = result || new Matrix();
    var m = result.m;

    m[0] = x;
    m[1] = 0;
    m[2] = 0;
    m[3] = 0;

    m[4] = 0;
    m[5] = y;
    m[6] = 0;
    m[7] = 0;

    m[8] = 0;
    m[9] = 0;
    m[10] = z;
    m[11] = 0;

    m[12] = 0;
    m[13] = 0;
    m[14] = 0;
    m[15] = 1;

    return result;
};

// ### GL.Matrix.translate(x, y, z[, result])
//
// This emulates the OpenGL function `glTranslate()`. You can optionally pass
// an existing matrix in `result` to avoid allocating a new matrix.
Matrix.translate = function (x, y, z, result) {
    result = result || new Matrix();
    var m = result.m;

    m[0] = 1;
    m[1] = 0;
    m[2] = 0;
    m[3] = x;

    m[4] = 0;
    m[5] = 1;
    m[6] = 0;
    m[7] = y;

    m[8] = 0;
    m[9] = 0;
    m[10] = 1;
    m[11] = z;

    m[12] = 0;
    m[13] = 0;
    m[14] = 0;
    m[15] = 1;

    return result;
};

// ### GL.Matrix.rotate(a, x, y, z[, result])
//
// Returns a matrix that rotates by `a` degrees around the vector `x, y, z`.
// You can optionally pass an existing matrix in `result` to avoid allocating
// a new matrix. This emulates the OpenGL function `glRotate()`.
Matrix.rotate = function (a, x, y, z, result) {
    if (!a || (!x && !y && !z)) {
        return Matrix.identity(result);
    }

    result = result || new Matrix();
    var m = result.m;

    var d = Math.sqrt(x * x + y * y + z * z);
    a *= Math.PI / 180; x /= d; y /= d; z /= d;
    var c = Math.cos(a), s = Math.sin(a), t = 1 - c;

    m[0] = x * x * t + c;
    m[1] = x * y * t - z * s;
    m[2] = x * z * t + y * s;
    m[3] = 0;

    m[4] = y * x * t + z * s;
    m[5] = y * y * t + c;
    m[6] = y * z * t - x * s;
    m[7] = 0;

    m[8] = z * x * t - y * s;
    m[9] = z * y * t + x * s;
    m[10] = z * z * t + c;
    m[11] = 0;

    m[12] = 0;
    m[13] = 0;
    m[14] = 0;
    m[15] = 1;

    return result;
};


//****************************************************************************************************************************************************************************************************************


// Represents indexed triangle geometry with arbitrary additional attributes.
// You need a shader to draw a mesh; meshes can't draw themselves.
//
// A mesh is a collection of `GL.Buffer` objects which are either vertex buffers
// (holding per-vertex attributes) or index buffers (holding the order in which
// vertices are rendered). By default, a mesh has a position vertex buffer called
// `vertices` and a triangle index buffer called `triangles`. New buffers can be
// added using `addVertexBuffer()` and `addIndexBuffer()`. Two strings are
// required when adding a new vertex buffer, the name of the data array on the
// mesh instance and the name of the GLSL attribute in the vertex shader.
//
// Example usage:
//
//     var mesh = new GL.Mesh({ coords: true, lines: true });
//
//     // Default attribute "vertices", available as "gl_Vertex" in
//     // the vertex shader
//     mesh.vertices = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0]];
//
//     // Optional attribute "coords" enabled in constructor,
//     // available as "gl_TexCoord" in the vertex shader
//     mesh.coords = [[0, 0], [1, 0], [0, 1], [1, 1]];
//
//     // Custom attribute "weights", available as "weight" in the
//     // vertex shader
//     mesh.addVertexBuffer('weights', 'weight');
//     mesh.weights = [1, 0, 0, 1];
//
//     // Default index buffer "triangles"
//     mesh.triangles = [[0, 1, 2], [2, 1, 3]];
//
//     // Optional index buffer "lines" enabled in constructor
//     mesh.lines = [[0, 1], [0, 2], [1, 3], [2, 3]];
//
//     // Upload provided data to GPU memory
//     mesh.compile();

// ### new GL.Indexer()
//
// Generates indices into a list of unique objects from a stream of objects
// that may contain duplicates. This is useful for generating compact indexed
// meshes from unindexed data.
function Indexer() {
    this.unique = [];
    this.indices = [];
    this.map = {};
}

Indexer.prototype = {
    // ### .add(v)
    //
    // Adds the object `obj` to `unique` if it hasn't already been added. Returns
    // the index of `obj` in `unique`.
    add: function (obj) {
        var key = JSON.stringify(obj);
        if (!(key in this.map)) {
            this.map[key] = this.unique.length;
            this.unique.push(obj);
        }
        return this.map[key];
    }
};

// ### new GL.Buffer(target, type)
//
// Provides a simple method of uploading data to a GPU buffer. Example usage:
//
//     var vertices = new GL.Buffer(gl.ARRAY_BUFFER, Float32Array);
//     var indices = new GL.Buffer(gl.ELEMENT_ARRAY_BUFFER, Uint16Array);
//     vertices.data = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0]];
//     indices.data = [[0, 1, 2], [2, 1, 3]];
//     vertices.compile();
//     indices.compile();
//
function Buffer(target, type) {
    this.buffer = null;
    this.target = target;
    this.type = type;
    this.data = [];
}

Buffer.prototype = {
    // ### .compile(type)
    //
    // Upload the contents of `data` to the GPU in preparation for rendering. The
    // data must be a list of lists where each inner list has the same length. For
    // example, each element of data for vertex normals would be a list of length three.
    // This will remember the data length and element length for later use by shaders.
    // The type can be either `gl.STATIC_DRAW` or `gl.DYNAMIC_DRAW`, and defaults to
    // `gl.STATIC_DRAW`.
    //
    // This could have used `[].concat.apply([], this.data)` to flatten
    // the array but Google Chrome has a maximum number of arguments so the
    // concatenations are chunked to avoid that limit.
    compile: function (type) {
        var data = [];
        for (var i = 0, chunk = 10000; i < this.data.length; i += chunk) {
            data = Array.prototype.concat.apply(data, this.data.slice(i, i + chunk));
        }
        var spacing = this.data.length ? data.length / this.data.length : 0;
        if (spacing !== Math.round(spacing)) throw new Error('not consistent size');
        this.buffer = this.buffer || gl.createBuffer();
        this.buffer.length = data.length;
        this.buffer.spacing = spacing;
        gl.bindBuffer(this.target, this.buffer);
        gl.bufferData(this.target, new this.type(data), type || gl.STATIC_DRAW);
    }
};

// ### new GL.Mesh([options])
//
// Represents a collection of vertex buffers and index buffers. Each vertex
// buffer maps to one attribute in GLSL and has a corresponding property set
// on the Mesh instance. There is one vertex buffer by default: `vertices`,
// which maps to `gl_Vertex`. The `coords`, `normals`, and `colors` vertex
// buffers map to `gl_TexCoord`, `gl_Normal`, and `gl_Color` respectively,
// and can be enabled by setting the corresponding options to true. There are
// two index buffers, `triangles` and `lines`, which are used for rendering
// `gl.TRIANGLES` and `gl.LINES`, respectively. Only `triangles` is enabled by
// default, although `computeWireframe()` will add a normal buffer if it wasn't
// initially enabled.
function Mesh(options) {
    options = options || {};
    this.vertexBuffers = {};
    this.indexBuffers = {};
    this.addVertexBuffer('vertices', 'gl_Vertex');
    this.dynamic = false;
    if (options.coords) this.addVertexBuffer('coords', 'gl_TexCoord');
    if (options.normals) this.addVertexBuffer('normals', 'gl_Normal');
    if (options.colors) this.addVertexBuffer('colors', 'gl_Color');
    if (!('triangles' in options) || options.triangles) this.addIndexBuffer('triangles');
    if (options.lines) this.addIndexBuffer('lines');
}

Mesh.prototype = {
    // ### .addVertexBuffer(name, attribute)
    //
    // Add a new vertex buffer with a list as a property called `name` on this object
    // and map it to the attribute called `attribute` in all shaders that draw this mesh.
    addVertexBuffer: function (name, attribute) {
        var buffer = this.vertexBuffers[attribute] = new Buffer(gl.ARRAY_BUFFER, Float32Array);
        buffer.name = name;
        this[name] = [];
    },

    // ### .addIndexBuffer(name)
    //
    // Add a new index buffer with a list as a property called `name` on this object.
    addIndexBuffer: function (name) {
        var buffer = this.indexBuffers[name] = new Buffer(gl.ELEMENT_ARRAY_BUFFER, Uint16Array);
        this[name] = [];
    },

    // ### .compile()
    //
    // Upload all attached buffers to the GPU in preparation for rendering. This
    // doesn't need to be called every frame, only needs to be done when the data
    // changes.
    compile: function () {
        for (var attribute in this.vertexBuffers) {
            var vertexBuffer = this.vertexBuffers[attribute];
            vertexBuffer.data = this[vertexBuffer.name];
            vertexBuffer.compile();
        }

        for (var name in this.indexBuffers) {
            var indexBuffer = this.indexBuffers[name];
            indexBuffer.data = this[name];
            indexBuffer.compile();
        }
    },

    // ### .transform(matrix)
    //
    // Transform all vertices by `matrix` and all normals by the inverse transpose
    // of `matrix`.
    transform: function (matrix) {
        var vertices = this.dynamic ? this.verticesOriginal : this.vertices;
        var normals = this.dynamic ? this.normalsOriginal : this.normals;

        this.vertices = vertices.map(function (v) {
            return matrix.transformPoint(Vector.fromArray(v)).toArray();
        });
        if (this.normals) {
            var invTrans = matrix.inverse().transpose();
            this.normals = normals.map(function (n) {
                return invTrans.transformVector(Vector.fromArray(n)).unit().toArray();
            });
        }
        this.compile();
        return this;
    },

    // ### .computeNormals()
    //
    // Computes a new normal for each vertex from the average normal of the
    // neighboring triangles. This means adjacent triangles must share vertices
    // for the resulting normals to be smooth.
    computeNormals: function () {
        if (!this.normals) this.addVertexBuffer('normals', 'gl_Normal');
        for (var i0 = 0; i0 < this.vertices.length; i0++) {
            this.normals[i0] = new Vector();
        }
        for (var i1 = 0; i1 < this.triangles.length; i1++) {
            var t = this.triangles[i1];
            var a = Vector.fromArray(this.vertices[t[0]]);
            var b = Vector.fromArray(this.vertices[t[1]]);
            var c = Vector.fromArray(this.vertices[t[2]]);
            var normal = b.subtract(a).cross(c.subtract(a)).unit();
            this.normals[t[0]] = this.normals[t[0]].add(normal);
            this.normals[t[1]] = this.normals[t[1]].add(normal);
            this.normals[t[2]] = this.normals[t[2]].add(normal);
        }
        for (var i2 = 0; i2 < this.vertices.length; i2++) {
            this.normals[i2] = this.normals[i2].unit().toArray();
        }
        this.compile();
        return this;
    },


    // ### .getAABB()
    //
    // Computes the axis-aligned bounding box, which is an object whose `min` and
    // `max` properties contain the minimum and maximum coordinates of all vertices.
    getAABB: function () {
        var aabb = { min: new Vector(Number.MAX_VALUE, Number.MAX_VALUE, Number.MAX_VALUE) };
        aabb.max = aabb.min.negative();
        for (var i = 0; i < this.vertices.length; i++) {
            var v = Vector.fromArray(this.vertices[i]);
            aabb.min = Vector.min(aabb.min, v);
            aabb.max = Vector.max(aabb.max, v);
        }
        return aabb;
    },

    createOriginalVerticesNormals: function () {
        this.dynamic = true;

        this.verticesOriginal = new Array();
        for (var indexVertices = 0; indexVertices < this.vertices.length; indexVertices++) {
            this.verticesOriginal.push([this.vertices[indexVertices][0], this.vertices[indexVertices][1], this.vertices[indexVertices][2]]);
        }

        this.normalsOriginal = new Array();
        for (var indexNormals = 0; indexNormals < this.normals.length; indexNormals++) {
            this.normalsOriginal.push([this.normals[indexNormals][0], this.normals[indexNormals][1], this.normals[indexNormals][2]]);
        }
    }
};

// ### GL.Mesh.plane([options])
//
// Generates a square 2x2 mesh the xy plane centered at the origin. The
// `options` argument specifies options to pass to the mesh constructor.
// Additional options include `detailX` and `detailY`, which set the tesselation
// in x and y, and `detail`, which sets both `detailX` and `detailY` at once.
// Two triangles are generated by default.
// Example usage:
//
//     var mesh1 = GL.Mesh.plane();
//     var mesh2 = GL.Mesh.plane({ detail: 5 });
//     var mesh3 = GL.Mesh.plane({ detailX: 20, detailY: 40 });
//
Mesh.plane = function (options) {
    options = options || {};
    var mesh = new Mesh(options);
    detailX = options.detailX || options.detail || 1;
    detailY = options.detailY || options.detail || 1;

    for (var y = 0; y <= detailY; y++) {
        var t = y / detailY;
        for (var x = 0; x <= detailX; x++) {
            var s = x / detailX;
            mesh.vertices.push([2 * s - 1, 2 * t - 1, 0]);
            if (mesh.coords) mesh.coords.push([s, t]);
            if (mesh.normals) mesh.normals.push([0, 0, 1]);
            if (x < detailX && y < detailY) {
                var i = x + y * (detailX + 1);
                mesh.triangles.push([i, i + 1, i + detailX + 1]);
                mesh.triangles.push([i + detailX + 1, i + 1, i + detailX + 2]);
            }
        }
    }

    mesh.compile();
    return mesh;
};

Mesh.square = function (p0, p1, p2, p3) {
    var mesh = new Mesh();
    mesh.vertices.push([p0.x, p0.y, p0.z]);
    mesh.vertices.push([p1.x, p1.y, p1.z]);
    mesh.vertices.push([p2.x, p2.y, p2.z]);
    mesh.vertices.push([p3.x, p3.y, p3.z]);
    mesh.triangles.push([0, 1, 2]);
    mesh.triangles.push([2, 1, 3]);
    mesh.compile();
    return mesh;
};

var cubeData = [
    [0, 4, 2, 6, -1, 0, 0], // -x
    [1, 3, 5, 7, +1, 0, 0], // +x
    [0, 1, 4, 5, 0, -1, 0], // -y
    [2, 6, 3, 7, 0, +1, 0], // +y
    [0, 2, 1, 3, 0, 0, -1], // -z
    [4, 5, 6, 7, 0, 0, +1]  // +z
];

function pickOctant(i) {
    return new Vector((i & 1) * 2 - 1, (i & 2) - 1, (i & 4) / 2 - 1);
}

// ### GL.Mesh.cube([options])
//
// Generates a 2x2x2 box centered at the origin. The `options` argument
// specifies options to pass to the mesh constructor.
Mesh.cube = function (options) {
    var mesh = new Mesh(options);

    for (var i = 0; i < cubeData.length; i++) {
        var data = cubeData[i], v = i * 4;
        for (var j = 0; j < 4; j++) {
            var d = data[j];
            mesh.vertices.push(pickOctant(d).toArray());
            if (mesh.coords) mesh.coords.push([j & 1, (j & 2) / 2]);
            if (mesh.normals) mesh.normals.push(data.slice(4, 7));
        }
        mesh.triangles.push([v, v + 1, v + 2]);
        mesh.triangles.push([v + 2, v + 1, v + 3]);
    }

    mesh.compile();
    return mesh;
};

// ### GL.Mesh.sphere([options])
//
// Generates a geodesic sphere of radius 1. The `options` argument specifies
// options to pass to the mesh constructor in addition to the `detail` option,
// which controls the tesselation level. The detail is `6` by default.
// Example usage:
//
//     var mesh1 = GL.Mesh.sphere();
//     var mesh2 = GL.Mesh.sphere({ detail: 2 });
//
Mesh.sphere = function (options) {
    function tri(a, b, c) { return flip ? [a, c, b] : [a, b, c]; }
    function fix(x) { return x + (x - x * x) / 2; }
    options = options || {};
    var mesh = new Mesh(options);
    var indexer = new Indexer();
    detail = options.detail || 6;

    for (var octant = 0; octant < 8; octant++) {
        var scale = pickOctant(octant);
        var flip = scale.x * scale.y * scale.z > 0;
        var data = [];
        for (var i = 0; i <= detail; i++) {
            // Generate a row of vertices on the surface of the sphere
            // using barycentric coordinates.
            for (var j0 = 0; i + j0 <= detail; j0++) {
                var a0 = i / detail;
                var b0 = j0 / detail;
                var c = (detail - i - j0) / detail;
                var vertex = { vertex: new Vector(fix(a0), fix(b0), fix(c)).unit().multiply(scale).toArray() };
                if (mesh.coords) vertex.coord = scale.y > 0 ? [1 - a0, c] : [c, 1 - a0];
                data.push(indexer.add(vertex));
            }

            // Generate triangles from this row and the previous row.
            if (i > 0) {
                for (var j1 = 0; i + j1 <= detail; j1++) {
                    var a1 = (i - 1) * (detail + 1) + ((i - 1) - (i - 1) * (i - 1)) / 2 + j1;
                    var b1 = i * (detail + 1) + (i - i * i) / 2 + j1;
                    mesh.triangles.push(tri(data[a1], data[a1 + 1], data[b1]));
                    if (i + j1 < detail) {
                        mesh.triangles.push(tri(data[b1], data[a1 + 1], data[b1 + 1]));
                    }
                }
            }
        }
    }

    // Reconstruct the geometry from the indexer.
    mesh.vertices = indexer.unique.map(function (v) { return v.vertex; });
    if (mesh.coords) mesh.coords = indexer.unique.map(function (v) { return v.coord; });
    if (mesh.normals) mesh.normals = mesh.vertices;

    mesh.compile();
    return mesh;
};



//****************************************************************************************************************************************************************************************************************


// Provides a convenient wrapper for WebGL shaders. A few uniforms and attributes,
// prefixed with `gl_`, are automatically added to all shader sources to make
// simple shaders easier to write.
//
// Example usage:
//
//     var shader = new GL.Shader('\
//       void main() {\
//         gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;\
//       }\
//     ', '\
//       uniform vec4 color;\
//       void main() {\
//         gl_FragColor = color;\
//       }\
//     ');
//
//     shader.uniforms({
//       color: [1, 0, 0, 1]
//     }).draw(mesh);

function regexMap(regex, text, callback) {
    while ((result = regex.exec(text)) !== null) {
        callback(result);
    }
}

// Non-standard names beginning with `gl_` must be mangled because they will
// otherwise cause a compiler error.
var LIGHTGL_PREFIX = 'LIGHTGL';

// ### new GL.Shader(vertexSource, fragmentSource)
//
// Compiles a shader program using the provided vertex and fragment shaders.
function Shader(vertexSource, fragmentSource) {
    // Allow passing in the id of an HTML script tag with the source
    function followScriptTagById(id) {
        var element = document.getElementById(id);
        return element ? element.text : id;
    }
    vertexSource = followScriptTagById(vertexSource);
    fragmentSource = followScriptTagById(fragmentSource);

    // Headers are prepended to the sources to provide some automatic functionality.
    var header = '\
    uniform mat3 gl_NormalMatrix;\
    uniform mat4 gl_ModelViewMatrix;\
    uniform mat4 gl_ProjectionMatrix;\
    uniform mat4 gl_ModelViewProjectionMatrix;\
    uniform mat4 gl_ModelViewMatrixInverse;\
    uniform mat4 gl_ProjectionMatrixInverse;\
    uniform mat4 gl_ModelViewProjectionMatrixInverse;\
  ';
    var vertexHeader = header + '\
    attribute vec4 gl_Vertex;\
    attribute vec4 gl_TexCoord;\
    attribute vec3 gl_Normal;\
    attribute vec4 gl_Color;\
    vec4 ftransform() {\
      return gl_ModelViewProjectionMatrix * gl_Vertex;\
    }\
  ';
    var fragmentHeader = '\
    precision highp float;\
  ' + header;

    // Check for the use of built-in matrices that require expensive matrix
    // multiplications to compute, and record these in `usedMatrices`.
    var source = vertexSource + fragmentSource;
    var usedMatrices = {};
    regexMap(/\b(gl_[^;]*)\b;/g, header, function (groups) {
        var name = groups[1];
        if (source.indexOf(name) !== -1) {
            var capitalLetters = name.replace(/[a-z_]/g, '');
            usedMatrices[capitalLetters] = LIGHTGL_PREFIX + name;
        }
    });
    if (source.indexOf('ftransform') !== -1) usedMatrices.MVPM = LIGHTGL_PREFIX + 'gl_ModelViewProjectionMatrix';
    this.usedMatrices = usedMatrices;

    // The `gl_` prefix must be substituted for something else to avoid compile
    // errors, since it's a reserved prefix. This prefixes all reserved names with
    // `_`. The header is inserted after any extensions, since those must come
    // first.
    function fix(header, source) {
        var replaced = {};
        var match = /^((\s*\/\/.*\n|\s*#extension.*\n)+)[^]*$/.exec(source);
        source = match ? match[1] + header + source.substr(match[1].length) : header + source;
        regexMap(/\bgl_\w+\b/g, header, function (result) {
            if (!(result in replaced)) {
                source = source.replace(new RegExp('\\b' + result + '\\b', 'g'), LIGHTGL_PREFIX + result);
                replaced[result] = true;
            }
        });
        return source;
    }
    vertexSource = fix(vertexHeader, vertexSource);
    fragmentSource = fix(fragmentHeader, fragmentSource);

    // Compile and link errors are thrown as strings.
    function compileSource(type, source) {
        var shader = gl.createShader(type);
        gl.shaderSource(shader, source);
        gl.compileShader(shader);
        if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS)) {
            throw new Error('compile error: ' + gl.getShaderInfoLog(shader));
        }
        return shader;
    }
    this.program = gl.createProgram();
    gl.attachShader(this.program, compileSource(gl.VERTEX_SHADER, vertexSource));
    gl.attachShader(this.program, compileSource(gl.FRAGMENT_SHADER, fragmentSource));
    gl.linkProgram(this.program);
    if (!gl.getProgramParameter(this.program, gl.LINK_STATUS)) {
        throw new Error('link error: ' + gl.getProgramInfoLog(this.program));
    }
    this.attributes = {};
    this.uniformLocations = {};

    // Sampler uniforms need to be uploaded using `gl.uniform1i()` instead of `gl.uniform1f()`.
    // To do this automatically, we detect and remember all uniform samplers in the source code.
    var isSampler = {};
    regexMap(/uniform\s+sampler(1D|2D|3D|Cube)\s+(\w+)\s*;/g, vertexSource + fragmentSource, function (groups) {
        isSampler[groups[2]] = 1;
    });
    this.isSampler = isSampler;
}

function isArray(obj) {
    var str = Object.prototype.toString.call(obj);
    return str === '[object Array]' || str === '[object Float32Array]';
}

function isNumber(obj) {
    var str = Object.prototype.toString.call(obj);
    return str === '[object Number]' || str === '[object Boolean]';
}

var tempMatrix = new Matrix();
var resultMatrix = new Matrix();

Shader.prototype = {
    // ### .uniforms(uniforms)
    //
    // Set a uniform for each property of `uniforms`. The correct `gl.uniform*()` method is
    // inferred from the value types and from the stored uniform sampler flags.
    uniforms: function (uniforms) {
        gl.useProgram(this.program);

        for (var name in uniforms) {
            var location = this.uniformLocations[name] || gl.getUniformLocation(this.program, name);
            if (!location) continue;
            this.uniformLocations[name] = location;
            var value = uniforms[name];
            if (value instanceof Vector) {
                value = [value.x, value.y, value.z];
            } else if (value instanceof Matrix) {
                value = value.m;
            }
            if (isArray(value)) {
                switch (value.length) {
                    case 1: gl.uniform1fv(location, new Float32Array(value)); break;
                    case 2: gl.uniform2fv(location, new Float32Array(value)); break;
                    case 3: gl.uniform3fv(location, new Float32Array(value)); break;
                    case 4: gl.uniform4fv(location, new Float32Array(value)); break;
                    // Matrices are automatically transposed, since WebGL uses column-major
                    // indices instead of row-major indices.
                    case 9: gl.uniformMatrix3fv(location, false, new Float32Array([
                        value[0], value[3], value[6],
                        value[1], value[4], value[7],
                        value[2], value[5], value[8]
                    ])); break;
                    case 16: gl.uniformMatrix4fv(location, false, new Float32Array([
                        value[0], value[4], value[8], value[12],
                        value[1], value[5], value[9], value[13],
                        value[2], value[6], value[10], value[14],
                        value[3], value[7], value[11], value[15]
                    ])); break;
                    default: throw new Error('don\'t know how to load uniform "' + name + '" of length ' + value.length);
                }
            } else if (isNumber(value)) {
                (this.isSampler[name] ? gl.uniform1i : gl.uniform1f).call(gl, location, value);
            } else {
                throw new Error('attempted to set uniform "' + name + '" to invalid value ' + value);
            }
        }

        return this;
    },

    // ### .draw(mesh[, mode])
    //
    // Sets all uniform matrix attributes, binds all relevant buffers, and draws the
    // mesh geometry as indexed triangles or indexed lines. Set `mode` to `gl.LINES`
    // (and either add indices to `lines` or call `computeWireframe()`) to draw the
    // mesh in wireframe.
    draw: function (mesh, mode) {
        this.drawBuffers(mesh.vertexBuffers,
            mesh.indexBuffers[mode === gl.LINES ? 'lines' : 'triangles'],
            arguments.length < 2 ? gl.TRIANGLES : mode);
    },

    // ### .drawBuffers(vertexBuffers, indexBuffer, mode)
    //
    // Sets all uniform matrix attributes, binds all relevant buffers, and draws the
    // indexed mesh geometry. The `vertexBuffers` argument is a map from attribute
    // names to `Buffer` objects of type `gl.ARRAY_BUFFER`, `indexBuffer` is a `Buffer`
    // object of type `gl.ELEMENT_ARRAY_BUFFER`, and `mode` is a WebGL primitive mode
    // like `gl.TRIANGLES` or `gl.LINES`. This method automatically creates and caches
    // vertex attribute pointers for attributes as needed.
    drawBuffers: function (vertexBuffers, indexBuffer, mode) {
        // Only construct up the built-in matrices we need for this shader.
        var used = this.usedMatrices;
        var MVM = gl.modelviewMatrix;
        var PM = gl.projectionMatrix;
        var MVMI = (used.MVMI || used.NM) ? MVM.inverse() : null;
        var PMI = (used.PMI) ? PM.inverse() : null;
        var MVPM = (used.MVPM || used.MVPMI) ? PM.multiply(MVM) : null;
        var matrices = {};
        if (used.MVM) matrices[used.MVM] = MVM;
        if (used.MVMI) matrices[used.MVMI] = MVMI;
        if (used.PM) matrices[used.PM] = PM;
        if (used.PMI) matrices[used.PMI] = PMI;
        if (used.MVPM) matrices[used.MVPM] = MVPM;
        if (used.MVPMI) matrices[used.MVPMI] = MVPM.inverse();
        if (used.NM) {
            var m = MVMI.m;
            matrices[used.NM] = [m[0], m[4], m[8], m[1], m[5], m[9], m[2], m[6], m[10]];
        }
        this.uniforms(matrices);

        // Create and enable attribute pointers as necessary.
        var length = 0;
        for (var vertex in vertexBuffers) {
            var buffer = vertexBuffers[vertex];
            var location = this.attributes[vertex] ||
                gl.getAttribLocation(this.program, vertex.replace(/^(gl_.*)$/, LIGHTGL_PREFIX + '$1'));
            if (location === -1 || !buffer.buffer) continue;
            this.attributes[vertex] = location;
            gl.bindBuffer(gl.ARRAY_BUFFER, buffer.buffer);
            gl.enableVertexAttribArray(location);
            gl.vertexAttribPointer(location, buffer.buffer.spacing, gl.FLOAT, false, 0, 0);
            length = buffer.buffer.length / buffer.buffer.spacing;
        }

        // Disable unused attribute pointers.
        for (var attribute in this.attributes) {
            if (!(attribute in vertexBuffers)) {
                gl.disableVertexAttribArray(this.attributes[attribute]);
            }
        }

        // Draw the geometry.
        if (length && (!indexBuffer || indexBuffer.buffer)) {
            if (indexBuffer) {
                gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, indexBuffer.buffer);
                gl.drawElements(mode, indexBuffer.buffer.length, gl.UNSIGNED_SHORT, 0);
            } else {
                gl.drawArrays(mode, 0, length);
            }
        }

        return this;
    }
};


//****************************************************************************************************************************************************************************************************************


// Provides a simple 3D vector class. Vector operations can be done using member
// functions, which return new vectors, or static functions, which reuse
// existing vectors to avoid generating garbage.
function Vector(x, y, z) {
    this.x = x || 0;
    this.y = y || 0;
    this.z = z || 0;
}

// ### Instance Methods
// The methods `add()`, `subtract()`, `multiply()`, and `divide()` can all
// take either a vector or a number as an argument.
Vector.prototype = {
    set: function (x, y, z) {
        this.x = x; this.y = y; this.z = z;
    },
    negative: function () {
        return new Vector(-this.x, -this.y, -this.z);
    },
    add: function (v) {
        if (v instanceof Vector) return new Vector(this.x + v.x, this.y + v.y, this.z + v.z);
        else return new Vector(this.x + v, this.y + v, this.z + v);
    },
    subtract: function (v) {
        if (v instanceof Vector) return new Vector(this.x - v.x, this.y - v.y, this.z - v.z);
        else return new Vector(this.x - v, this.y - v, this.z - v);
    },
    multiply: function (v) {
        if (v instanceof Vector) return new Vector(this.x * v.x, this.y * v.y, this.z * v.z);
        else return new Vector(this.x * v, this.y * v, this.z * v);
    },
    divide: function (v) {
        if (v instanceof Vector) return new Vector(this.x / v.x, this.y / v.y, this.z / v.z);
        else return new Vector(this.x / v, this.y / v, this.z / v);
    },
    equals: function (v) {
        return this.x === v.x && this.y === v.y && this.z === v.z;
    },
    dot: function (v) {
        return this.x * v.x + this.y * v.y + this.z * v.z;
    },
    cross: function (v) {
        return new Vector(
            this.y * v.z - this.z * v.y,
            this.z * v.x - this.x * v.z,
            this.x * v.y - this.y * v.x
        );
    },
    length: function () {
        return Math.sqrt(this.dot(this));
    },
    unit: function () {
        return this.divide(this.length());
    },
    min: function () {
        return Math.min(Math.min(this.x, this.y), this.z);
    },
    max: function () {
        return Math.max(Math.max(this.x, this.y), this.z);
    },
    toArray: function (n) {
        return [this.x, this.y, this.z].slice(0, n || 3);
    },
    clone: function () {
        return new Vector(this.x, this.y, this.z);
    },
    init: function (x, y, z) {
        this.x = x; this.y = y; this.z = z;
        return this;
    }
};

// ### Static Methods
// `Vector.randomDirection()` returns a vector with a length of 1 and a
// statistically uniform direction. `Vector.lerp()` performs linear
// interpolation between two vectors.
Vector.negative = function (a, b) {
    b.x = -a.x; b.y = -a.y; b.z = -a.z;
    return b;
};
Vector.add = function (a, b, c) {
    if (b instanceof Vector) { c.x = a.x + b.x; c.y = a.y + b.y; c.z = a.z + b.z; }
    else { c.x = a.x + b; c.y = a.y + b; c.z = a.z + b; }
    return c;
};
Vector.subtract = function (a, b, c) {
    if (b instanceof Vector) { c.x = a.x - b.x; c.y = a.y - b.y; c.z = a.z - b.z; }
    else { c.x = a.x - b; c.y = a.y - b; c.z = a.z - b; }
    return c;
};
Vector.scale = function (a, b) {
    a.x *= b;
    a.y *= b;
    a.z *= b;
};
Vector.multiply = function (a, b, c) {
    if (b instanceof Vector) { c.x = a.x * b.x; c.y = a.y * b.y; c.z = a.z * b.z; }
    else { c.x = a.x * b; c.y = a.y * b; c.z = a.z * b; }
    return c;
};
Vector.divide = function (a, b, c) {
    if (b instanceof Vector) { c.x = a.x / b.x; c.y = a.y / b.y; c.z = a.z / b.z; }
    else { c.x = a.x / b; c.y = a.y / b; c.z = a.z / b; }
    return c;
};
Vector.cross = function (a, b, c) {
    c.x = a.y * b.z - a.z * b.y;
    c.y = a.z * b.x - a.x * b.z;
    c.z = a.x * b.y - a.y * b.x;
    return c;
};
Vector.unit = function (a, b) {
    var length = a.length();
    b.x = a.x / length;
    b.y = a.y / length;
    b.z = a.z / length;
    return b;
};
Vector.fromAngles = function (theta, phi) {
    return new Vector(Math.cos(theta) * Math.cos(phi), Math.sin(phi), Math.sin(theta) * Math.cos(phi));
};
Vector.min = function (a, b) {
    return new Vector(Math.min(a.x, b.x), Math.min(a.y, b.y), Math.min(a.z, b.z));
};
Vector.max = function (a, b) {
    return new Vector(Math.max(a.x, b.x), Math.max(a.y, b.y), Math.max(a.z, b.z));
};
Vector.lerp = function (a, b, fraction) {
    if (fraction > 1.0)
        return b;
    return b.subtract(a).multiply(fraction).add(a);
};
Vector.fromArray = function (a) {
    return new Vector(a[0], a[1], a[2]);
};
Vector.distance = function (v1, v2) {
    var dx = v1.x - v2.x;
    var dy = v1.y - v2.y;
    var dz = v1.z - v2.z;

    return Math.sqrt(dx * dx + dy * dy + dz * dz);
}


//****************************************************************************************************************************************************************************************************************



// The internal `gl` variable holds the current WebGL context.
var gl;

var GL = {
    // ### Initialization
    //
    // `GL.create()` creates a new WebGL context and augments it with more
    // methods. The alpha channel is disabled by default because it usually causes
    // unintended transparencies in the canvas.
    create: function (options) {
        options = options || {};
        var canvas = document.createElement('canvas');
        canvas.width = 800;
        canvas.height = 600;
        if (!('alpha' in options)) options.alpha = false;
        try { gl = canvas.getContext('webgl', options); } catch (ex) { console.log(ex) }
        try { gl = canvas.getContext('experimental-webgl', options); } catch (ex) { console.log(ex) }
        if (!gl) throw new Error('WebGL not supported');
        gl.HALF_FLOAT_OES = 0x8D61;
        addMatrixStack();
        addEventListeners();
        addOtherMethods();
        return gl;
    },

    // `GL.keys` contains a mapping of key codes to booleans indicating whether
    // that key is currently pressed.
    keys: {},

    // Export all external classes.
    Matrix: Matrix,
    Indexer: Indexer,
    Buffer: Buffer,
    Mesh: Mesh,
    //HitTest: HitTest,
    //Raytracer: Raytracer,
    Shader: Shader,
    //Texture: Texture,
    Vector: Vector
};

// ### Matrix stack
//
// Implement the OpenGL modelview and projection matrix stacks, along with some
// other useful GLU matrix functions.

function addMatrixStack() {
    gl.MODELVIEW = ENUM | 1;
    gl.PROJECTION = ENUM | 2;
    var tempMatrix = new Matrix();
    var resultMatrix = new Matrix();
    gl.modelviewMatrix = new Matrix();
    gl.projectionMatrix = new Matrix();
    var modelviewStack = [];
    var projectionStack = [];
    var matrix, stack;
    gl.matrixMode = function (mode) {
        switch (mode) {
            case gl.MODELVIEW:
                matrix = 'modelviewMatrix';
                stack = modelviewStack;
                break;
            case gl.PROJECTION:
                matrix = 'projectionMatrix';
                stack = projectionStack;
                break;
            default:
                throw new Error('invalid matrix mode ' + mode);
        }
    };
    gl.loadIdentity = function () {
        Matrix.identity(gl[matrix]);
    };
    gl.loadMatrix = function (m) {
        var from = m.m, to = gl[matrix].m;
        for (var i = 0; i < 16; i++) {
            to[i] = from[i];
        }
    };
    gl.multMatrix = function (m) {
        gl.loadMatrix(Matrix.multiply(gl[matrix], m, resultMatrix));
    };
    gl.perspective = function (fov, aspect, near, far) {
        gl.multMatrix(Matrix.perspective(fov, aspect, near, far, tempMatrix));
    };
    gl.frustum = function (l, r, b, t, n, f) {
        gl.multMatrix(Matrix.frustum(l, r, b, t, n, f, tempMatrix));
    };
    gl.scale = function (x, y, z) {
        gl.multMatrix(Matrix.scale(x, y, z, tempMatrix));
    };
    gl.translate = function (x, y, z) {
        gl.multMatrix(Matrix.translate(x, y, z, tempMatrix));
    };
    gl.rotate = function (a, x, y, z) {
        gl.multMatrix(Matrix.rotate(a, x, y, z, tempMatrix));
    };
    //gl.lookAt = function (ex, ey, ez, cx, cy, cz, ux, uy, uz) {
    //    gl.multMatrix(Matrix.lookAt(ex, ey, ez, cx, cy, cz, ux, uy, uz, tempMatrix));
    //};
    //gl.pushMatrix = function () {
    //    stack.push(Array.prototype.slice.call(gl[matrix].m));
    //};
    //gl.popMatrix = function () {
    //    var m = stack.pop();
    //    gl[matrix].m = hasFloat32Array ? new Float32Array(m) : m;
    //};
    //gl.project = function (objX, objY, objZ, modelview, projection, viewport) {
    //    modelview = modelview || gl.modelviewMatrix;
    //    projection = projection || gl.projectionMatrix;
    //    viewport = viewport || gl.getParameter(gl.VIEWPORT);
    //    var point = projection.transformPoint(modelview.transformPoint(new Vector(objX, objY, objZ)));
    //    return new Vector(
    //        viewport[0] + viewport[2] * (point.x * 0.5 + 0.5),
    //        viewport[1] + viewport[3] * (point.y * 0.5 + 0.5),
    //        point.z * 0.5 + 0.5
    //    );
    //};
    //gl.unProject = function (winX, winY, winZ, modelview, projection, viewport) {
    //    modelview = modelview || gl.modelviewMatrix;
    //    projection = projection || gl.projectionMatrix;
    //    viewport = viewport || gl.getParameter(gl.VIEWPORT);
    //    var point = new Vector(
    //        (winX - viewport[0]) / viewport[2] * 2 - 1,
    //        (winY - viewport[1]) / viewport[3] * 2 - 1,
    //        winZ * 2 - 1
    //    );
    //    return Matrix.inverse(Matrix.multiply(projection, modelview, tempMatrix), resultMatrix).transformPoint(point);
    //};
    gl.matrixMode(gl.MODELVIEW);
}

// ### Improved mouse events
//
// This adds event listeners on the `gl.canvas` element that call
// `gl.onmousedown()`, `gl.onmousemove()`, and `gl.onmouseup()` with an
// augmented event object. The event object also has the properties `x`, `y`,
// `deltaX`, `deltaY`, and `dragging`.
function addEventListeners() {
    var context = gl, oldX = 0, oldY = 0, buttons = {}, hasOld = false;
    var has = Object.prototype.hasOwnProperty;
    function isDragging() {
        if (isLocked)
            return true;
        for (var b in buttons) {
            if (has.call(buttons, b) && buttons[b]) return true;
        }
        return false;
    }
    function augment(original) {
        // Make a copy of original, a native `MouseEvent`, so we can overwrite
        // WebKit's non-standard read-only `x` and `y` properties (which are just
        // duplicates of `pageX` and `pageY`). We can't just use
        // `Object.create(original)` because some `MouseEvent` functions must be
        // called in the context of the original event object.
        var e = {};
        for (var name in original) {
            if (typeof original[name] === 'function') {
                e[name] = (function (callback) {
                    return function () {
                        callback.apply(original, arguments);
                    };
                })(original[name]);
            } else {
                e[name] = original[name];
            }
        }
        e.original = original;
        e.x = e.pageX;
        e.y = e.pageY;
        for (var obj = gl.canvas; obj; obj = obj.offsetParent) {
            e.x -= obj.offsetLeft;
            e.y -= obj.offsetTop;
        }
        if (hasOld) {
            e.deltaX = e.x - oldX;
            e.deltaY = e.y - oldY;
        } else {
            e.deltaX = 0;
            e.deltaY = 0;
            hasOld = true;
        }
        oldX = e.x;
        oldY = e.y;
        e.dragging = isDragging();
        if (e.dragging) {
            e.deltaX = e.movementX;
            e.deltaY = e.movementY;
        }
        e.preventDefault = function () {
            e.original.preventDefault();
        };
        e.stopPropagation = function () {
            e.original.stopPropagation();
        };

        return e;
    }
    function mousedown(e) {
        gl = context;
        if (!isDragging()) {
            // Expand the event handlers to the document to handle dragging off canvas.
            on(document, 'mousemove', mousemove);
            on(document, 'mouseup', mouseup);
            off(gl.canvas, 'mousemove', mousemove);
            off(gl.canvas, 'mouseup', mouseup);
        }
        buttons[e.which] = true;
        e = augment(e);
        if (gl.onmousedown) gl.onmousedown(e);
        e.preventDefault();
    }
    function mousemove(e) {
        gl = context;
        e = augment(e);
        if (gl.onmousemove) gl.onmousemove(e);
        e.preventDefault();
    }
    function mouseup(e) {
        gl = context;
        buttons[e.which] = false;
        if (!isDragging()) {
            // Shrink the event handlers back to the canvas when dragging ends.
            off(document, 'mousemove', mousemove);
            off(document, 'mouseup', mouseup);
            on(gl.canvas, 'mousemove', mousemove);
            on(gl.canvas, 'mouseup', mouseup);
        }
        e = augment(e);
        if (gl.onmouseup) gl.onmouseup(e);
        e.preventDefault();
    }
    function reset() {
        hasOld = false;
    }
    function resetAll() {
        buttons = {};
        hasOld = false;
    }
    on(gl.canvas, 'mousedown', mousedown);
    on(gl.canvas, 'mousemove', mousemove);
    on(gl.canvas, 'mouseup', mouseup);
    on(gl.canvas, 'mouseover', reset);
    on(gl.canvas, 'mouseout', reset);
    on(document, 'contextmenu', resetAll);
}

// ### Automatic keyboard state
//
// The current keyboard state is stored in `GL.keys`, a map of integer key
// codes to booleans indicating whether that key is currently pressed. Certain
// keys also have named identifiers that can be used directly, such as
// `GL.keys.SPACE`. Values in `GL.keys` are initially undefined until that
// key is pressed for the first time. If you need a boolean value, you can
// cast the value to boolean by applying the not operator twice (as in
// `!!GL.keys.SPACE`).

function mapKeyCode(code) {
    var named = {
        8: 'BACKSPACE',
        9: 'TAB',
        13: 'ENTER',
        16: 'SHIFT',
        27: 'ESCAPE',
        32: 'SPACE',
        37: 'LEFT',
        38: 'UP',
        39: 'RIGHT',
        40: 'DOWN'
    };
    return named[code] || (code >= 65 && code <= 90 ? String.fromCharCode(code) : null);
}

function on(element, name, callback) {
    element.addEventListener(name, callback);
}

function off(element, name, callback) {
    element.removeEventListener(name, callback);
}

on(document, 'keydown', function (e) {
    if (!e.altKey && !e.ctrlKey && !e.metaKey) {
        var key = mapKeyCode(e.keyCode);
        if (key) GL.keys[key] = true;
        GL.keys[e.keyCode] = true;
    }
});

on(document, 'keyup', function (e) {
    if (!e.altKey && !e.ctrlKey && !e.metaKey) {
        var key = mapKeyCode(e.keyCode);
        if (key) GL.keys[key] = false;
        GL.keys[e.keyCode] = false;
    }
});

function addOtherMethods() {
    // ### Multiple contexts
    //
    // When using multiple contexts in one web page, `gl.makeCurrent()` must be
    // called before issuing commands to a different context.
    (function (context) {
        gl.makeCurrent = function () {
            gl = context;
        };
    })(gl);

    // ### Animation
    //
    // Call `gl.animate()` to provide an animation loop that repeatedly calls
    // `gl.onupdate()` and `gl.ondraw()`.
    gl.animate = function () {
        var time = new Date().getTime();
        var context = gl;

        function update() {
            gl = context;
            var now = new Date().getTime();
            if (gl.onupdate) gl.onupdate((now - time) / 1000);
            if (gl.ondraw) gl.ondraw();
            time = now;
            requestAnimationFrame(update);
        }

        update();
    };

    // ### Fullscreen
    //
    // Provide an easy way to get a fullscreen app running, including an
    // automatic 3D perspective projection matrix by default. This should be
    // called once.
    //
    // Just fullscreen, no automatic camera:
    //
    //     gl.fullscreen({ camera: false });
    //
    // Adjusting field of view, near plane distance, and far plane distance:
    //
    //     gl.fullscreen({ fov: 45, near: 0.1, far: 1000 });
    //
    // Adding padding from the edge of the window:
    //
    //     gl.fullscreen({ paddingLeft: 250, paddingBottom: 60 });
    //
    gl.fullscreen = function (options) {
        options = options || {};
        var top = options.paddingTop || 0;
        var left = options.paddingLeft || 0;
        var right = options.paddingRight || 0;
        var bottom = options.paddingBottom || 0;
        if (!document.body) {
            throw new Error('gl.fullscreen()');
        }
        document.body.appendChild(gl.canvas);
        document.body.style.overflow = 'hidden';
        gl.canvas.style.position = 'absolute';
        gl.canvas.style.left = left + 'px';
        gl.canvas.style.top = top + 'px';
        function resize() {
            gl.canvas.width = window.innerWidth - left - right;
            gl.canvas.height = window.innerHeight - top - bottom;
            gl.viewport(0, 0, gl.canvas.width, gl.canvas.height);
            if (options.camera || !('camera' in options)) {
                gl.matrixMode(gl.PROJECTION);
                gl.loadIdentity();
                gl.perspective(options.fov || 45, gl.canvas.width / gl.canvas.height,
                    options.near || 0.1, options.far || 1000);
                gl.matrixMode(gl.MODELVIEW);
            }
            if (gl.ondraw) gl.ondraw();
        }
        on(window, 'resize', resize);
        resize();
    };
}

// A value to bitwise-or with new enums to make them distinguishable from the
// standard WebGL enums.
var ENUM = 0x12340000;

function getRandom(min, max) {
    return Math.random() * (max - min) + min;
}

function getRandomInt(min, max) {
    min = Math.ceil(min);
    max = Math.floor(max);
    return Math.floor(Math.random() * (max - min)) + min; //The maximum is exclusive and the minimum is inclusive
}

function degToRad(degrees) {
    return degrees * Math.PI / 180;
}

function intersection(ro, rd, aabb) {
    var out = new Float32Array(3);
    var d = distance(ro, rd, aabb);
    if (d === +Infinity || d === -Infinity) {
        out = null;
    } else {
        out = out || [];
        for (var i = 0; i < ro.length; i++) {
            out[i] = ro[i] + rd[i] * d;
        }
    }
    if (out === null)
        return null;
    else
        return new GL.Vector(out[0], out[1], out[2]);
}

function distance(ro, rd, aabb) {
    var dims = ro.length;
    var lo = -Infinity;
    var hi = +Infinity;

    for (var i = 0; i < dims; i++) {
        var p1 = 0; var p2 = 0;
        if (i === 0) {
            p1 = aabb.min.x;
            p2 = aabb.max.x;
        }
        if (i === 1) {
            p1 = aabb.min.y;
            p2 = aabb.max.y;
        }
        if (i === 2) {
            p1 = aabb.min.z;
            p2 = aabb.max.z;
        }
        var dimLo = (p1 - ro[i]) / rd[i];
        var dimHi = (p2 - ro[i]) / rd[i];

        if (dimLo > dimHi) {
            var tmp = dimLo;
            dimLo = dimHi;
            dimHi = tmp;
        }

        if (dimHi < lo || dimLo > hi) {
            return Infinity;
        }

        if (dimLo > lo) lo = dimLo;
        if (dimHi < hi) hi = dimHi;
    }

    return lo > hi ? Infinity : lo;
}

function inteceptCircleLineSeg(circle, p1, p2, radius) {
    var a, b, c, d, u1, u2, ret, retP1, retP2, v1, v2;
    v1 = {};
    v2 = {};
    v1.x = p2.z - p1.z;
    v1.y = p2.x - p1.x;
    v2.x = p1.z - circle.z;
    v2.y = p1.x - circle.x;
    b = (v1.x * v2.x + v1.y * v2.y);
    c = 2 * (v1.x * v1.x + v1.y * v1.y);
    b *= -2;
    d = Math.sqrt(b * b - 2 * c * (v2.x * v2.x + v2.y * v2.y - radius * radius));
    if (isNaN(d)) { // no intercept
        return false;
    }
    u1 = (b - d) / c;  // these represent the unit distance of point one and two on the line
    u2 = (b + d) / c;
    retP1 = {};   // return points
    retP2 = {}
    ret = []; // return array
    if (u1 <= 1 && u1 >= 0) {  // add point if on the line segment
        retP1.x = p1.z + v1.x * u1;
        retP1.y = p1.x + v1.y * u1;
        ret[0] = retP1;
    }
    if (u2 <= 1 && u2 >= 0) {  // second add point if on the line segment
        retP2.x = p1.z + v1.x * u2;
        retP2.y = p1.x + v1.y * u2;
        ret[ret.length] = retP2;
    }
    if (ret.length === 2 || ret.length === 1)
        return true;
    else
        return false;
}

function newMaze(x, y) {
    // Establish variables and starting grid
    var totalCells = x * y;
    var cells = new Array();
    var unvis = new Array();
    for (var i = 0; i < y; i++) {
        cells[i] = new Array();
        unvis[i] = new Array();
        for (var j = 0; j < x; j++) {
            cells[i][j] = [0, 0, 0, 0];
            unvis[i][j] = true;
        }
    }

    // Set a random position to start from
    var currentCell = [Math.floor(Math.random() * y), Math.floor(Math.random() * x)];
    var path = [currentCell];
    unvis[currentCell[0]][currentCell[1]] = false;
    var visited = 1;

    // Loop through all available cell positions
    while (visited < totalCells) {
        // Determine neighboring cells
        var pot = [[currentCell[0] - 1, currentCell[1], 0, 2],
        [currentCell[0], currentCell[1] + 1, 1, 3],
        [currentCell[0] + 1, currentCell[1], 2, 0],
        [currentCell[0], currentCell[1] - 1, 3, 1]];
        var neighbors = new Array();

        // Determine if each neighboring cell is in game grid, and whether it has already been checked
        for (var l = 0; l < 4; l++) {
            if (pot[l][0] > -1 && pot[l][0] < y && pot[l][1] > -1 && pot[l][1] < x && unvis[pot[l][0]][pot[l][1]]) { neighbors.push(pot[l]); }
        }

        // If at least one active neighboring cell has been found
        if (neighbors.length) {
            // Choose one of the neighbors at random
            next = neighbors[Math.floor(Math.random() * neighbors.length)];

            // Remove the wall between the current cell and the chosen neighboring cell
            cells[currentCell[0]][currentCell[1]][next[2]] = 1;
            cells[next[0]][next[1]][next[3]] = 1;

            // Mark the neighbor as visited, and set it as the current cell
            unvis[next[0]][next[1]] = false;
            visited++;
            currentCell = [next[0], next[1]];
            path.push(currentCell);
        }
        // Otherwise go back up a step and keep going
        else {
            currentCell = path.pop();
        }
    }
    return cells;
}

function pointInAABB(point, aabb) {
    //Check if the point is less than max and greater than min
    if (point.x > aabb.min.x && point.x < aabb.max.x &&
        point.y > aabb.min.y && point.y < aabb.max.y &&
        point.z > aabb.min.z && point.z < aabb.max.z) {
        return true;
    }

    //If not, then return false
    return false;
}

function lerp(value1, value2, amount) {
    amount = amount < 0 ? 0 : amount;
    amount = amount > 1 ? 1 : amount;
    return value1 + (value2 - value1) * amount;
}

function easeInSine(timeElapsed, startValue, c, totalTime) {
    return -c * Math.cos(timeElapsed / totalTime * (Math.PI / 2)) + c + startValue;
}

function removeIsDeadArray(array) {

    if (array === null)
        return;

    if (array.length === 0)
        return;

    var indexesToRemove = new Array();

    for (var indexArray = 0; indexArray < array.length; indexArray++)
        if (array[indexArray].isDead())
            indexesToRemove.push(indexArray);

    for (var index = 0; index < indexesToRemove.length; index++)
        array.splice(indexesToRemove[index], 1);
}

function getCookie(cname) {
    var name = cname + "=";
    var decodedCookie = decodeURIComponent(document.cookie);
    var ca = decodedCookie.split(';');
    for (var i = 0; i < ca.length; i++) {
        var c = ca[i];
        while (c.charAt(0) === ' ') {
            c = c.substring(1);
        }
        if (c.indexOf(name) === 0) {
            return c.substring(name.length, c.length);
        }
    }
    return "";
}


//********************************************************************************************************************************************************************************************************************************


function Hero(arrayBoxEnemies) {
    this.Camera = new GL.Vector(0.0, 0.0, 0.0);
    this.arrayEffects = new Array();
    this.arrayBoxEnemies = arrayBoxEnemies;
    this.health = 1.0;

    this.update = function (step) {

        for (var indexEffects = 0; indexEffects < this.arrayEffects.length; indexEffects++)
            this.arrayEffects[indexEffects].update(step);

        removeIsDeadArray(this.arrayEffects);

        var minDistToEnemyScanning = Number.MAX_SAFE_INTEGER;
        for (var indexEnemies = 0; indexEnemies < this.arrayBoxEnemies.length; indexEnemies++) {
            if (this.arrayBoxEnemies[indexEnemies].isDead() === false && this.arrayBoxEnemies[indexEnemies].isScanActive) {
                minDistToEnemyScanning = Math.min(minDistToEnemyScanning, GL.Vector.distance(this.Camera, this.arrayBoxEnemies[indexEnemies].position));
            }
        }
        if (minDistToEnemyScanning < 10.0)
            playScanSoundEffect(minDistToEnemyScanning);

    }
    this.playScanSoundEffect = function (distance) {
    }
    this.takeHit = function () {        
        this.arrayEffects.push(new HeroCameraShakeEfect(this, 0.5, 12000));       
        this.health -= 0.20;
    }
    this.isDead = function () {
        if (this.health < 0.0) return true;
        else return false;
    }
    this.addCameraShakeEffect = function () {
        this.arrayEffects.push(new HeroCameraShakeEfect(this, 0.15, 6000));
    }
}


function HeroCameraShakeEfect(hero, lifespan, stepMultiplier) {
    this.hero = hero;
    this.lifespan = lifespan;
    this.joggingAngle = 0;

    this.update = function (step) {
        this.lifespan -= step;
        this.joggingAngle += step * stepMultiplier;
        this.hero.Camera.y = 2.0 + Math.sin(degToRad(this.joggingAngle)) / 20;
    }
    this.isDead = function () {
        if (this.lifespan < 0) return true;
        else return false;
    }
}


//****************************************************************************************************************************************************************************************************************

function BoxEnemy(startPosition, arrayMaze, blockSize, oneBlockSpeed, scanTotalTime, scanRadius, hero) {
    this.arrayMaze = arrayMaze;
    this.position = startPosition;
    this.lastDirInverted = -1;
    this.effects = new Array();
    this.fireSpheres = new Array();
    this.blockSize = blockSize;
    this.oneBlockSpeed = oneBlockSpeed;
    this.scanTotalTime = scanTotalTime;
    this.scanRadius = scanRadius;
    this.isScanActive = false;
    this.health = 1.0;
    this.color = new GL.Vector(0.0, this.health, 0.0);
    this.hero = hero;
    this.mesh = GL.Mesh.cube({ normals: true, coords: true  }).transform(GL.Matrix.scale(1.0, 1.0, 1.0)).transform(GL.Matrix.translate(this.position.x, this.position.y, this.position.z));
    //this.mesh.createOriginalVerticesNormals();


    this.update = function (step) {
        if (this.effects.length === 0) {
            this.effects.push(new ScanEffect(this, this.scanTotalTime));
        }


        for (var indexEffects = 0; indexEffects < this.effects.length; indexEffects++) {
            this.effects[indexEffects].update(step);
        }
        removeIsDeadArray(this.effects);


        for (var indexFireSpheres = 0; indexFireSpheres < this.fireSpheres.length; indexFireSpheres++) {
            this.fireSpheres[indexFireSpheres].update(step);
        }
        removeIsDeadArray(this.fireSpheres);

        //this.mesh = this.mesh.transform(GL.Matrix.translate(this.position.x, this.position.y, this.position.z));
        this.mesh = GL.Mesh.cube({ normals: true, coords: true }).transform(GL.Matrix.scale(1.0, 1.0, 1.0)).transform(GL.Matrix.translate(this.position.x, this.position.y, this.position.z));

        var distanceToHero = GL.Vector.distance(this.hero.Camera, this.position);
        if (distanceToHero < 2.0) {
            this.hero.takeHit();
        }
    }

    this.isDead = function () {
        if (this.health <= 0)
            return true;
        else
            return false;
    }

    this.isHeroNear = function () {
        var distance = GL.Vector.distance(this.position, this.hero.Camera);
        if (distance < this.scanRadius) return true;
        else return false;
    }

    this.fire = function () {
        this.fireSpheres.push(new FireSpheres(this, this.scanRadius, this.scanTotalTime));
    }

    this.takeHit = function () {
        this.health -= 0.25;
        this.color = this.color.subtract(0.25);
    }

    this.nextNormalMovement = function () {
        var currentMazeBlock = null;

        for (var indexObjects = 0; indexObjects < arrayMaze.length; indexObjects++) {
            if (pointInAABB(this.position, arrayMaze[indexObjects].AABB)) {
                currentMazeBlock = arrayMaze[indexObjects];
                break;
            }
        }

        if (currentMazeBlock !== null) {
            var finalDir = -1;

            var possibleDir = new Array();
            if (currentMazeBlock.MazePosition[0] === 1)
                possibleDir.push(0);
            if (currentMazeBlock.MazePosition[1] === 1)
                possibleDir.push(1);
            if (currentMazeBlock.MazePosition[2] === 1)
                possibleDir.push(2);
            if (currentMazeBlock.MazePosition[3] === 1)
                possibleDir.push(3);

            if (possibleDir.length === 1)
                finalDir = possibleDir[0];
            else {
                finalDir = possibleDir[getRandomInt(0, possibleDir.length)];
                for (var i = 0; i < 5; i++) {
                    if (this.lastDirInverted === finalDir)
                        finalDir = possibleDir[getRandomInt(0, possibleDir.length)];
                    else
                        break;
                }
            }

            if (finalDir !== -1) {
                var endPosition = null;
                if (finalDir === 0) {
                    endPosition = new GL.Vector(this.position.x - this.blockSize, this.position.y, this.position.z);
                    this.lastDirInverted = 2;
                }
                if (finalDir === 1) {
                    endPosition = new GL.Vector(this.position.x, this.position.y, this.position.z - this.blockSize);
                    this.lastDirInverted = 3;
                }
                if (finalDir === 2) {
                    endPosition = new GL.Vector(this.position.x + this.blockSize, this.position.y, this.position.z);
                    this.lastDirInverted = 0;
                }
                if (finalDir === 3) {
                    endPosition = new GL.Vector(this.position.x, this.position.y, this.position.z + this.blockSize);
                    this.lastDirInverted = 1;
                }

                if (endPosition !== null)
                    this.effects.push(new BoxMovementEfect(this, this.position, endPosition, this.oneBlockSpeed));
            }
        }
    }
}

function BoxMovementEfect(box, startPostion, endPosition, totalTime) {
    this.startPosition = startPostion.clone();
    this.endPosition = endPosition.clone();
    this.totalTime = totalTime;
    this.timeElapsed = 0.0;
    this.finished = false;

    this.update = function (step) {
        this.timeElapsed += step;
        box.position = GL.Vector.lerp(this.startPosition, this.endPosition, this.timeElapsed / this.totalTime);
    }

    this.isDead = function () {
        if (this.timeElapsed > this.totalTime) return true;
        else return false;
    }
}

function ScanEffect(box, totalTime) {
    this.box = box;
    this.totalTime = totalTime;
    this.timeElapsed = 0.0;
    this.boxOldColor = box.color.clone();
    this.box.color = new GL.Vector(0.0, 0.0, 0.0);
    this.up = true;
    this.timeElapsedColor = 0.0;

    this.update = function (step) {
        this.timeElapsed += step;
        this.timeElapsedColor += step;
        if (this.timeElapsed > this.totalTime) {

            if (this.box.isHeroNear()) {
                this.box.fire();
                this.timeElapsed = 0.0;
            }
            else {
                this.box.color = this.boxOldColor;
                this.box.nextNormalMovement();
            }
        }
        else {
            if (this.up === true)
                this.box.color.x = lerp(0.0, 1.0, this.timeElapsedColor);
            else
                this.box.color.x = lerp(1.0, 0.0, this.timeElapsedColor);

            if (this.timeElapsedColor > 1.0) {
                this.up = !this.up;
                this.timeElapsedColor = 0;
            }
        }

    }

    this.isDead = function () {
        if (this.timeElapsed > this.totalTime) return true;
        else return false;
    }
}


function FireSpheres(box, radius, totalTime) {
    this.box = box;
    this.totalTime = totalTime;
    this.timeElapsed = 0.0;
    this.maxRadius = radius;
    this.currentRadius = 0.0;
    this.mesh = GL.Mesh.sphere({ normals: true, coords: true }).transform(GL.Matrix.scale(this.currentRadius, this.currentRadius, this.currentRadius)).transform(GL.Matrix.translate(this.box.position.x, this.box.position.y, this.box.position.z));
    //this.mesh.createOriginalVerticesNormals();

    this.update = function (step) {
        this.timeElapsed += step;
        this.currentRadius = easeInSine(this.timeElapsed, 0.0, this.maxRadius, this.totalTime);
        this.mesh = GL.Mesh.sphere({ normals: true, coords: true }).transform(GL.Matrix.scale(this.currentRadius, this.currentRadius, this.currentRadius)).transform(GL.Matrix.translate(this.box.position.x, this.box.position.y, this.box.position.z));

        var distanceToHero = GL.Vector.distance(this.box.hero.Camera, this.box.position);
        if (distanceToHero < this.currentRadius) {
            this.box.hero.takeHit();
            this.timeElapsed = this.totalTime + step;
        }
    }

    this.isDead = function () {
        if (this.timeElapsed > this.totalTime)
            return true;
        else
            return false;
    }
}


//****************************************************************************************************************************************************************************************************************


function HitLine(p0, p1, p2, p3, lifespan) {
    this.mesh = new Mesh.square(p0, p1, p2, p3);
    this.lifespan = lifespan;

    this.update = function (step) {
        this.lifespan -= step;
    }

    this.isDead = function () {
        if (this.lifespan < 0) return true;
        else return false;
    }
}


//****************************************************************************************************************************************************************************************************************

function MazeBlock() {
    this.Center = null;
    this.AABB = null;
    this.MazePosition = null;

    this.ArrayObjects = new Array();
    this.ArrayAABBs = new Array();

    this.init = function (mazePosition, center) {
        this.MazePosition = mazePosition;
        this.Center = new GL.Vector(center.x, center.y, center.z * -1);

        var size = 2.0;
        var physicsSize = 1.95;

        var groundMesh = GL.Mesh.plane({ normals: true, coords: true }).transform(GL.Matrix.scale(size, size)).transform(GL.Matrix.rotate(-90.0, 1.0, 0.0, 0.0)).transform(GL.Matrix.translate(this.Center.x, 0.0, this.Center.z)).computeNormals();
        this.ArrayObjects.push(groundMesh);       
        this.ArrayAABBs.push(groundMesh.getAABB());

        if (mazePosition[0] === 0) {
            var xMinusMesh = GL.Mesh.plane({ normals: true, coords: true }).transform(GL.Matrix.scale(size, size)).transform(GL.Matrix.rotate(90.0, 0.0, 1.0, 0.0)).transform(GL.Matrix.translate(this.Center.x - size, 2.0, this.Center.z)).computeNormals();
            this.ArrayObjects.push(xMinusMesh);
            var meshXMinusMesh = GL.Mesh.plane({ normals: true, coords: true }).transform(GL.Matrix.scale(physicsSize, physicsSize)).transform(GL.Matrix.rotate(90.0, 0.0, 1.0, 0.0)).transform(GL.Matrix.translate(this.Center.x - physicsSize, 2.0, this.Center.z)).computeNormals();
            var xMinusAABB = meshXMinusMesh.getAABB();
            xMinusAABB.Plane = 0;
            this.ArrayAABBs.push(xMinusAABB);
        }

        if (mazePosition[1] === 0) {
            var zMinusMesh = GL.Mesh.plane({ normals: true, coords: true }).transform(GL.Matrix.scale(size, size)).transform(GL.Matrix.translate(this.Center.x, 2.0, (this.Center.z - size))).computeNormals();
            this.ArrayObjects.push(zMinusMesh);
            var meshZMinusMesh = GL.Mesh.plane({ normals: true, coords: true }).transform(GL.Matrix.scale(physicsSize, physicsSize)).transform(GL.Matrix.translate(this.Center.x, 2.0, (this.Center.z - physicsSize))).computeNormals();
            var zMinusAABB = meshZMinusMesh.getAABB();
            zMinusAABB.Plane = 1;
            this.ArrayAABBs.push(zMinusAABB);
        }

        if (mazePosition[2] === 0) {
            var xPlusMesh = GL.Mesh.plane({ normals: true, coords: true }).transform(GL.Matrix.scale(size, size)).transform(GL.Matrix.rotate(-90.0, 0.0, 1.0, 0.0)).transform(GL.Matrix.translate(this.Center.x + size, 2.0, this.Center.z)).computeNormals();
            this.ArrayObjects.push(xPlusMesh);
            var meshXPlusMesh = GL.Mesh.plane({ normals: true, coords: true }).transform(GL.Matrix.scale(physicsSize, physicsSize)).transform(GL.Matrix.rotate(-90.0, 0.0, 1.0, 0.0)).transform(GL.Matrix.translate(this.Center.x + physicsSize, 2.0, this.Center.z)).computeNormals();
            var xPlusAABB = meshXPlusMesh.getAABB();
            xPlusAABB.Plane = 2;
            this.ArrayAABBs.push(xPlusAABB);
        }

        if (mazePosition[3] === 0) {
            var zPlusMesh = GL.Mesh.plane({ normals: true, coords: true }).transform(GL.Matrix.scale(size, size)).transform(GL.Matrix.rotate(180.0, 0.0, 1.0, 0.0)).transform(GL.Matrix.translate(this.Center.x, 2.0, this.Center.z + size)).computeNormals();
            this.ArrayObjects.push(zPlusMesh);
            var meshZPlusMesh = GL.Mesh.plane({ normals: true, coords: true }).transform(GL.Matrix.scale(physicsSize, physicsSize)).transform(GL.Matrix.rotate(180.0, 0.0, 1.0, 0.0)).transform(GL.Matrix.translate(this.Center.x, 2.0, this.Center.z + physicsSize)).computeNormals();
            var zPlusAABB = meshZPlusMesh.getAABB();
            zPlusAABB.Plane = 3;
            this.ArrayAABBs.push(zPlusAABB);
        }

        var meshCube = GL.Mesh.cube().transform(GL.Matrix.scale(size, size, size)).transform(GL.Matrix.translate(this.Center.x, 2.0, this.Center.z));
        this.AABB = meshCube.getAABB();
    };
}


//****************************************************************************************************************************************************************************************************************


function Particle(options) {
    this.loc = options.loc || new GL.Vector();
    this.vel = options.vec || new GL.Vector();
    this.acc = options.acc || new GL.Vector();
    this.lifespan = options.lifespan || Math.random() * 5;

    this.update = function (step) {
        this.vel = this.vel.add(this.acc.multiply(step));
        this.loc = this.loc.add(this.vel.multiply(step));
        this.lifespan -= step;
    }

    this.isDead = function () {
        if (this.lifespan < 0) return true;
        else return false;
    }
}

//****************************************************************************************************************************************************************************************************************

var isLocked = false;

var arrayBullets = null;
var arrayParticles = null;
var arrayBoxEnemies = null;
var hitLineLifeSpan = 0.15;
var arrayHitLines = null;
var blockSize = 4.0;

var arrayMaze = null;
var fristMazeBlock = null;
var lastMazeBlock = null;

var currentLevel = 1;
var currentLevelIsReady = false;

var angleX = 0.0;
var angleY = 0.0;
var joggingAngle = 0;
var hero = null;
var gunCoolDown = 1.0;
var meshParticles = null;
var lastStep = 0.0;


var gl = GL.create();

// Event listener when the pointerlock is updated (or removed by pressing ESC for example).
var pointerlockchange = function () {
    var controlEnabled = document.mozPointerLockElement || document.webkitPointerLockElement || document.msPointerLockElement || document.pointerLockElement || null;
    // If the user is already locked
    if (!controlEnabled)
        isLocked = false;
    else
        isLocked = true;
};


var createParticlesAndLaserForShot = function (hitPoint) {
    var mat = new GL.Matrix();
    mat = GL.Matrix.rotate(angleY, 0, 1, 0, mat);
    var p0 = mat.transformVector(new GL.Vector(-0.05, -0.5, 0.0));
    var p1 = mat.transformVector(new GL.Vector(0.05, -0.5, 0.0));

    arrayHitLines.push(new HitLine(hero.Camera.add(p0), hero.Camera.add(p1), hitPoint, hitPoint, hitLineLifeSpan));
    hero.addCameraShakeEffect();

    for (var i = 0; i < 10; i++) {
        var particle = new Particle({ acc: new GL.Vector(0.0, -9.81, 0.0), vec: new GL.Vector(getRandom(-2.0, 2.0), getRandom(-2.0, 2.0), getRandom(-2.0, 2.0)), loc: hitPoint, lifespan: 1.5 });
        arrayParticles.push(particle);
    }
}

var loadCurrentLevel = function () {

    currentLevelIsReady = false;

    arrayBullets = new Array();
    arrayParticles = new Array();
    arrayBoxEnemies = new Array();
    arrayHitLines = new Array();

    arrayMaze = new Array();
    fristMazeBlock = null;
    lastMazeBlock = null;
    
    angleX = 0.0;
    angleY = 0.0;
    joggingAngle = 0;

    hero = new Hero(arrayBoxEnemies);
    hero.Camera = new GL.Vector(0.0, 2.0, 0.0);

    mazeCL = 3 + this.currentLevel;

    mazeObj = newMaze(mazeCL, mazeCL);

    for (var y = 0; y < mazeCL; y++) {
        for (var x = 0; x < mazeCL; x++) {
            var mazeBlock = new MazeBlock();
            mazeBlock.init(mazeObj[y][x], new GL.Vector(y * blockSize, -1.0, x * blockSize));
            arrayMaze.push(mazeBlock);
        }
    }

    fristMazeBlock = arrayMaze[0];
    lastMazeBlock = arrayMaze[arrayMaze.length - 1];

    var numEnemies = (mazeCL / 3) >> 0;
    for (var i = 0; i < numEnemies; i++) {
        var mazeBlockIndex = getRandomInt(1, arrayMaze.length);
        var speed = getRandom(2.0, 8.0);
        var scanTotalTime = getRandom(2.0, 8.0);
        var scanRadius = getRandom(5.0, 10.0);
        arrayBoxEnemies.push(new BoxEnemy(new GL.Vector(arrayMaze[mazeBlockIndex].Center.x, 2.0, arrayMaze[mazeBlockIndex].Center.z), arrayMaze, blockSize, speed, scanTotalTime, scanRadius, hero));
    }

    meshPortal = GL.Mesh.sphere({ normals: true, coords: true }).transform(GL.Matrix.scale(2.00, 2.00, 2.00)).transform(GL.Matrix.translate(lastMazeBlock.Center.x, 2.0, lastMazeBlock.Center.z)).computeNormals();
    meshPortalAABB = meshPortal.getAABB();

    var levelDescription = document.getElementById("levelDescription");
    levelDescription.style.display = "block";
    if (currentLevel === 1) {        
        levelDescription.innerHTML = "Welcome to Lost in the Maze<br/>Navigate the Maze and find the portal to the next Maze until you are free<br/>Be careful because the evil boxes can scan you and fire at you<br/>Good Luck";
    }
    else {
        levelDescription.innerHTML = "Level " + currentLevel;
    }
    document.cookie = "currentLevel=" + currentLevel;

    setTimeout(function () { var levelDescription = document.getElementById("levelDescription"); levelDescription.style.display = "none"; }, 5000);

    currentLevelIsReady = true;
}


function onloadFunction() {

    // Attach events to the document
    document.addEventListener("pointerlockchange", pointerlockchange, false);
    document.addEventListener("mspointerlockchange", pointerlockchange, false);
    document.addEventListener("mozpointerlockchange", pointerlockchange, false);
    document.addEventListener("webkitpointerlockchange", pointerlockchange, false);

    meshParticles = GL.Mesh.cube({ normals: true, coords: true }).transform(GL.Matrix.scale(0.05, 0.05, 0.05));
    meshParticles.createOriginalVerticesNormals();

    var currentLevelCookie = getCookie("currentLevel");
    if (currentLevelCookie !== null && currentLevelCookie !== "") {
        var currentLevelParsed = parseInt(currentLevelCookie);
        if (!isNaN(currentLevelParsed))
            currentLevel = currentLevelParsed;
    }   

    loadCurrentLevel();
}


var shader = new GL.Shader('vertex-Maze', 'fragment-Maze');
var shaderEnemies = new GL.Shader('vertex-Enemy', 'fragment-Enemy');
var shaderLines = new GL.Shader('vertex-id', 'fragment-id');
var shaderPortal = new GL.Shader('vertex-Portal', 'fragment-Portal');
var shaderBullets = new GL.Shader('vertex-Bullets', 'fragment-Bullets');

gl.canvas.onclick = function () {
    if (!isLocked) {
        gl.canvas.requestPointerLock = gl.canvas.requestPointerLock || gl.canvas.mozRequestPointerLock || gl.canvas.webkitRequestPointerLock;
        gl.canvas.requestPointerLock();
    }
    else {

        if (gunCoolDown > 0.0)
            return;

        var forward = GL.Vector.fromAngles((90 - angleY) * Math.PI / 180, (180 - angleX) * Math.PI / 180);
        var cameraForward = hero.Camera.add(forward.multiply(2.0));
        var direction = new GL.Vector(0, 0, 0);
        GL.Vector.subtract(hero.Camera, cameraForward, direction);

        const origin = new Float32Array([hero.Camera.x, hero.Camera.y, hero.Camera.z]);
        const dir = new Float32Array([direction.x, direction.y, direction.z]);

        var closestDist = Number.MAX_SAFE_INTEGER;
        var closestHitPoint = null;
        var closestHitPointEnemy = null;
        var enemyHit = null;
        for (var indexObjects = 0; indexObjects < arrayMaze.length; indexObjects++) {
            var rtBox = intersection(origin, dir, arrayMaze[indexObjects].AABB);
            if (rtBox !== null) {
                for (var indexAABBs = 0; indexAABBs < arrayMaze[indexObjects].ArrayAABBs.length; indexAABBs++) {
                    var hitPoint = intersection(origin, dir, arrayMaze[indexObjects].ArrayAABBs[indexAABBs]);
                    if (hitPoint !== null) {
                        var vectorCameraForward = hero.Camera.subtract(cameraForward);
                        var vectorCameraHitPoint = hero.Camera.subtract(hitPoint);
                        var dot = vectorCameraForward.dot(vectorCameraHitPoint);
                        //if (isNaN(angle) || (!isNaN(angle) && angle < (Math.PI * 0.5))) {
                        if (dot > 0.0) {
                            var dist = GL.Vector.distance(hero.Camera, hitPoint);
                            if (dist < closestDist) {
                                closestDist = dist;
                                closestHitPoint = hitPoint;
                            }
                        }
                    }
                }
            }
        }

        for (var indexEnemies = 0; indexEnemies < arrayBoxEnemies.length; indexEnemies++) {
            var hitPointEnemy = intersection(origin, dir, arrayBoxEnemies[indexEnemies].mesh.getAABB());
            if (hitPointEnemy !== null) {
                var vectorCameraForwardEnemy = hero.Camera.subtract(cameraForward);
                var vectorCameraHitPointEnemy = hero.Camera.subtract(hitPointEnemy);
                var dotEnemy = vectorCameraForwardEnemy.dot(vectorCameraHitPointEnemy);
                if (dotEnemy > 0.0) {
                    var distEnemy = GL.Vector.distance(hero.Camera, hitPointEnemy);
                    if (distEnemy < closestDist) {
                        closestDist = distEnemy;
                        closestHitPoint = null;
                        closestHitPointEnemy = hitPointEnemy;
                        enemyHit = arrayBoxEnemies[indexEnemies];
                    }
                }
            }
        }       

        if (closestHitPoint !== null) {
            var bullet = GL.Mesh.sphere({ normals: true, coords: true }).transform(GL.Matrix.scale(0.25, 0.25, 0.25)).transform(GL.Matrix.translate(closestHitPoint.x, closestHitPoint.y, closestHitPoint.z)).computeNormals();
            arrayBullets.push(bullet);
            createParticlesAndLaserForShot(closestHitPoint);
            gunCoolDown = 1.0;
        }
        else if (enemyHit !== null && closestHitPointEnemy !== null) {
            enemyHit.takeHit();
            createParticlesAndLaserForShot(closestHitPointEnemy);
            gunCoolDown = 1.0;
        }
    }
}

gl.onmousemove = function (e) {
    if (isLocked || e.dragging) {
        angleY -= e.deltaX * 0.25;
        angleX = Math.max(-90, Math.min(90, angleX - e.deltaY * 0.25));
    }
};

gl.onupdate = function (step) {

    lastStep += step;

    if (!currentLevelIsReady)
        return;

    if (gunCoolDown > 0.0)
        gunCoolDown -= step;

    var speed = step * 6;
    var newCamera = hero.Camera.clone();

    // Forward movement
    var up = GL.keys.W | GL.keys.UP;
    var down = GL.keys.S | GL.keys.DOWN;
    var forward = GL.Vector.fromAngles((90 - angleY) * Math.PI / 180, (180 - angleX) * Math.PI / 180);
    newCamera = newCamera.add(forward.multiply(speed * (up - down)));

    // Sideways movement
    var left = GL.keys.A | GL.keys.LEFT;
    var right = GL.keys.D | GL.keys.RIGHT;
    var sideways = GL.Vector.fromAngles(-angleY * Math.PI / 180, 0);
    newCamera = newCamera.add(sideways.multiply(speed * (right - left)));

    //Particles
    for (var indexParticles = 0; indexParticles < arrayParticles.length; indexParticles++) {
        arrayParticles[indexParticles].update(step);
    }
    removeIsDeadArray(arrayParticles);

    for (var indexHitLines = 0; indexHitLines < arrayHitLines.length; indexHitLines++) {
        arrayHitLines[indexHitLines].update(step);
    }
    removeIsDeadArray(arrayHitLines);

    //Enemies
    for (var indexEnemies = 0; indexEnemies < arrayBoxEnemies.length; indexEnemies++) {
        arrayBoxEnemies[indexEnemies].update(step);
    }
    removeIsDeadArray(arrayBoxEnemies);

    //Check for WIN
    if (pointInAABB(hero.Camera, meshPortalAABB)) {
        currentLevel++;
        loadCurrentLevel();
        return;
    }

    //Check GameOver
    if (hero.isDead()) {
        currentLevelIsReady = false;
        var levelDescription = document.getElementById("levelDescription");
        levelDescription.style.display = "block";        
        levelDescription.innerHTML = "HASTED";
        setTimeout(function () { var levelDescription = document.getElementById("levelDescription"); levelDescription.style.display = "none"; loadCurrentLevel(); }, 5000);
        return;
    }

    var hitPlaneZ = false;
    var hitPlaneX = false;

    for (var indexObjects = 0; indexObjects < arrayMaze.length; indexObjects++) {
        if (GL.Vector.distance(newCamera, arrayMaze[indexObjects].Center) <= 10.0) {
            for (var indexAABBs = 0; indexAABBs < arrayMaze[indexObjects].ArrayAABBs.length; indexAABBs++) {
                if (inteceptCircleLineSeg(newCamera, arrayMaze[indexObjects].ArrayAABBs[indexAABBs].min, arrayMaze[indexObjects].ArrayAABBs[indexAABBs].max, 1.0)) {
                    if (arrayMaze[indexObjects].ArrayAABBs[indexAABBs].Plane === 0)
                        hitPlaneX = true;
                    if (arrayMaze[indexObjects].ArrayAABBs[indexAABBs].Plane === 1)
                        hitPlaneZ = true;
                    if (arrayMaze[indexObjects].ArrayAABBs[indexAABBs].Plane === 2)
                        hitPlaneX = true;
                    if (arrayMaze[indexObjects].ArrayAABBs[indexAABBs].Plane === 3)
                        hitPlaneZ = true;
                }
            }
        }
    }

    if (hitPlaneZ)
        newCamera.z = hero.Camera.z;

    if (hitPlaneX)
        newCamera.x = hero.Camera.x;

    hero.Camera = newCamera;

    if (up || down || left || right) {
        joggingAngle += step * 600.0;
        hero.Camera.y = 2.0 + Math.sin(degToRad(joggingAngle)) / 20;
    }

    hero.update(step);
};

gl.ondraw = function () {

    if (!currentLevelIsReady)
        return;

    gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
    gl.loadIdentity();
    gl.rotate(-angleX, 1, 0, 0);
    gl.rotate(-angleY, 0, 1, 0);
    gl.translate(-hero.Camera.x, -hero.Camera.y, -hero.Camera.z);


    var distanceTotal = GL.Vector.distance(fristMazeBlock.Center, lastMazeBlock.Center);
    var distance = GL.Vector.distance(hero.Camera, lastMazeBlock.Center);

    var brig = 100.0 - ((100.0 * distance) / distanceTotal);
    if (brig < 10.0)
        brig = 10.0;

    //Maze
    for (var indexObjects = 0; indexObjects < arrayMaze.length; indexObjects++) {
        for (var indexMeshs = 0; indexMeshs < arrayMaze[indexObjects].ArrayObjects.length; indexMeshs++)
            shader.uniforms({ brightness: brig }).draw(arrayMaze[indexObjects].ArrayObjects[indexMeshs], gl.TRIANGLES);
    }

    //Bullets
    for (var indexBullets = 0; indexBullets < arrayBullets.length; indexBullets++) {
        shaderBullets.uniforms({ time: lastStep }).draw(arrayBullets[indexBullets], gl.TRIANGLES);
    }

    //Particles
    for (var indexParticles = 0; indexParticles < arrayParticles.length; indexParticles++) {
        meshParticles.transform(GL.Matrix.translate(arrayParticles[indexParticles].loc.x, arrayParticles[indexParticles].loc.y, arrayParticles[indexParticles].loc.z));
        shaderEnemies.uniforms({ time: lastStep, color: new GL.Vector(1.0, 1.0, 1.0) }).draw(meshParticles, gl.TRIANGLES);
    }

    //Enemies
    for (var indexEnemies = 0; indexEnemies < arrayBoxEnemies.length; indexEnemies++) {
        shaderEnemies.uniforms({ time: lastStep, color: arrayBoxEnemies[indexEnemies].color }).draw(arrayBoxEnemies[indexEnemies].mesh, gl.TRIANGLES);
    }

    //Enemies Spheres
    for (var indexEnemiesFireSpheres = 0; indexEnemiesFireSpheres < arrayBoxEnemies.length; indexEnemiesFireSpheres++) {
        for (var indexFireSpheres = 0; indexFireSpheres < arrayBoxEnemies[indexEnemiesFireSpheres].fireSpheres.length; indexFireSpheres++)
            shaderEnemies.uniforms({ time: lastStep, color: arrayBoxEnemies[indexEnemiesFireSpheres].color }).draw(arrayBoxEnemies[indexEnemiesFireSpheres].fireSpheres[indexFireSpheres].mesh, gl.TRIANGLES);
    }

    //HitLines
    for (var indexHitLines = 0; indexHitLines < arrayHitLines.length; indexHitLines++) {
        shaderBullets.uniforms({ brightness: brig }).draw(arrayHitLines[indexHitLines].mesh, gl.TRIANGLES);
    }

    //EndLevel Mesh
    shaderPortal.uniforms({ time: lastStep }).draw(meshPortal, gl.TRIANGLES);
};

gl.fullscreen();
gl.animate();
gl.enable(gl.CULL_FACE);
gl.enable(gl.POLYGON_OFFSET_FILL);
gl.polygonOffset(1, 1);
gl.clearColor(0.8, 0.8, 0.8, 1);
gl.enable(gl.DEPTH_TEST);

