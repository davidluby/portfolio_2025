'use client'
import React, { useEffect } from 'react'

const Sim = () => {
  const U_FIELD = 0
  const V_FIELD = 1
  const S_FIELD = 2
  let _count = 0

  const scene = {
    gravity : 0,
    dt : 1 / 60,
    iterations : 10,
    frame_nr : 0,
    over_relaxation : 1,
    obstacle_x : 0,
    obstacle_y : 0,
    obstacle_r : .05
  }

  useEffect (() => {
    // init canvas/context
    const canvas = document.getElementById('fluid')
    const ctx = canvas.getContext('2d')

    // device width/height is aspect ratio and limit drawn pixels
    const aspect = window.innerWidth / window.innerHeight
    canvas.width = 1000 * aspect
    canvas.height = 1000 / aspect


    // set vertical resolution and generate horizontal using aspect
    const resolution = 150
    const h = 1 / resolution

    const x_cells = Math.floor(aspect / h)
    const y_cells = Math.floor(1 / h)

    let flu = new Fluid(1000, x_cells, y_cells, h)

    update()

    function update () {
      if (scene.frame_nr % 100 == 0) {
        const x = (0.5 + Math.cos(Math.random() * Math.PI) / 4) * aspect
        const y = 0.5 + Math.cos(Math.random() * Math.PI) / 4
        interact(x, y, false)
      }
      flu.simulate(scene.dt, scene.gravity, scene.iterations)
      scene.frame_nr++
      draw()
      requestAnimationFrame(update)
    }

    function interact (x, y, reset) {
      console.log(x, y)
      let vx = 0
      let vy = 0

      if (!reset) {
        vx = (x - scene.obstacle_x) / scene.dt
        vy = (y - scene.obstacle_y) / scene.dt
      }

      scene.obstacle_x = x
      scene.obstacle_y = y
      const r = scene.obstacle_r
      const n = flu.y_dim

      for (let i = 1; i < flu.x_dim - 2; i++) {
        for (let j = 1; j < flu.y_dim - 2; j++) {
          flu.s[i*n +j] = 1

          let dx = (i + 0.5) * flu.h - x
          let dy = (j + 0.5) * flu.h - y

          if (dx * dx + dy * dy < r * r) {
            flu.s[i*n + j] = 1.0
            flu.d[i*n + j] = 0.5 + Math.sin(0.1 * scene.frame_nr) / 2
            flu.u[i*n + j] = vx
            flu.u[(i+1)*n + j] = vx
            flu.v[i*n + j] = vy
            flu.v[i*n + j+1] = vy
          }
        }
      }
    }

    function draw () {
      // clear canvas initialize image data (array of pixels)
      ctx.clearRect(0, 0, canvas.width, canvas.height)
      let id = ctx.getImageData(0, 0, canvas.width, canvas.height)
  
      // length and width of drawn fluid element
      const x_step = canvas.width / (flu.x_dim - 4)
      const y_step = canvas.height / (flu.y_dim - 4)


      // dont draw boundary cells (for 0 < j < y_dim)
      for (let i = 2; i < flu.x_dim - 2; i++) {
        for (let j = 2; j < flu.y_dim - 2; j++) {
          let col = color(flu.d[i * flu.y_dim + j])

          // offset x and y by -1 due to boundary cells
          let x = Math.floor((i - 2) * x_step)
          let y = Math.floor((j - 2) * y_step)

          // for height of element (in pixels)
          for (let yi = y; yi < y + y_step; yi++) {
            // go to the zeroth pixel column and the yith pixel row
            let p = 4 * (yi * canvas.width + x)

            // for width of element (in pixels)
            for (let xi = 0; xi < x_step; xi++) {
              // go to the xith pixel column at the yith pixel row and update rgba
              id.data[p++] = col[0]
              id.data[p++] = col[1]
              id.data[p++] = col[2]
              id.data[p++] = col[3]
            }
          }
        }
      }
      ctx.putImageData(id, 0, 0)
    }
    
    function color (val) {
      /* const min = 0
      const max = 1

      val = Math.min(Math.max(val, min), max - 0.0001)
      const diff = max - min
      val = diff == 0.0 ? 0.5 : (val - min) / diff

      const m = 0.25
      const num = Math.floor(val / m)
      const s = (val - num * m) / m
      let r, b, g
      switch (num) {
        case 0 : r = 0.0; g = s; b = 1; break
        case 1 : r = 0.0; g = 1.0; b = 1.0 - s; break;
        case 2 : r = s; g = 1.0; b = 0.0; break;
        case 3 : r = 1.0; g = 1.0 - s; b = 0.0; break
      } */

      //return [val * 255, val * 255, val * 255, 255]
      return [val * 21, val * 76, val * 121, 255]
      //return [r * 255, g * 255, b * 255, 255]
    }

    var mouseDown = false;

	function startDrag(x, y) {
		mouseDown = true;
		x = x / window.innerWidth * aspect
    y = y / window.innerHeight
		interact(x, y, true);
	}

	function drag(x, y) {
		if (mouseDown) {
			x = x / window.innerWidth * aspect
      y = y / window.innerHeight
			interact(x, y, false);
		}
	}

	function endDrag() {
		mouseDown = false;
	}

	canvas.addEventListener('mousedown', event => {
		startDrag(event.x, event.y);
	});

	canvas.addEventListener('mouseup', _event => {
		endDrag();
	});

	canvas.addEventListener('mousemove', event => {
		drag(event.x, event.y);
	});

	canvas.addEventListener('touchstart', event => {
		startDrag(event.touches[0].clientX, event.touches[0].clientY)
	});

	canvas.addEventListener('touchend', _event => {
		endDrag()
	});
  canvas.addEventListener('touchmove', event => {
		event.preventDefault();
		event.stopImmediatePropagation();
		drag(event.touches[0].clientX, event.touches[0].clientY)
	}, { passive: false});
  }, [])

    class Fluid {
      constructor (rho, x_dim, y_dim, h) {
        // fluid density, x/y resolutions, element height as percentage
        this.rho = rho
        this.x_dim = x_dim + 2
        this.y_dim = y_dim + 2
        this.h = h

        // fluid "density" or cell fullness, pressure/smoke fields
        this.d = new Float32Array(this.x_dim * this.y_dim)
        this.d.fill(1.0)
        this.d_new = new Float32Array(this.x_dim * this.y_dim)
        this.p = new Float32Array(this.x_dim * this.y_dim)
        this.s = new Float32Array(this.x_dim * this.y_dim)

        // horizontal/vertical velocity right/up positive
        this.u = new Float32Array(this.x_dim * this.y_dim)
        this.v = new Float32Array(this.x_dim * this.y_dim)
        this.u_new = new Float32Array(this.x_dim * this.y_dim)
        this.v_new = new Float32Array(this.x_dim * this.y_dim)

        //
      }


      // simulation routine
      integrate (dt, gravity) {
        const n = this.y_dim
        for (let i = 1; i < this.x_dim; i++ ) {
          for (let j = 1; j < this.y_dim - 1; j++) {
            if (this.s[i*n + j] != 0.0 && this.s[i*n + j - 1] != 0.0) {
              this.v[i*n + j] += gravity * dt
            }
          }
        }
      }

      solve_incompressibility (iterations, dt) {
        const n = this.y_dim
        const cp = this.rho * this.h / dt

        for (let k = 0; k < iterations; k++) {
          for (let i = 1; i < this.x_dim - 1; i++) {
            for (let j = 1; j < this.y_dim - 1; j++) {
              
              if (this.s[i*n + j] == 0.0) {
                continue
              }
              // in = out
              let s = this.s[i*n + j]
              let sx0 = this.s[(i-1)*n + j]
              let sx1 = this.s[(i+1)*n + j]
              let sy0 = this.s[i*n + j-1]
              let sy1 = this.s[i*n + j+1]
              s = sx0 + sx1 + sy0 + sy1
              // confirm
              if (s == 0.0) {
                continue
              }

              let div = this.u[(i+1)*n + j] - this.u[i*n + j] +
                        this.v[i*n + j+1] - this.v[i*n + j]

              let p = -div / s
              p *= scene.over_relaxation
              this.p[i*n + j] += cp * p

              this.u[i*n + j] -= sx0 * p
              this.u[(i+1)*n + j] += sx1 * p
              this.v[i*n + j] -= sy0 * p
              this.v[i*n + j+1] += sy1 * p              
            }
          }
        }
      }

      extrapolate () {
        const n = this.y_dim
        for (let i = 0; i < this.x_dim; i++) {
          this.u[i*n + 0] = this.u[i*n + 1]
          this.u[i*n + this.y_dim - 1] = this.u[i*n + this.y_dim - 2]
        }

        for (let j = 0; j < this.y_dim; j++) {
          this.v[0*n + j] = this.v[1*n + j]
          this.v[(this.x_dim-1)*n + j] = this.v[(this.x_dim-2)*n + j]
        }
      }

      sample_field (x, y, field) {
        const n = this.y_dim
        const h1 = 1 / this.h
        const h2 = this.h / 2

        x = Math.max(Math.min(x, this.x_dim * this.h), this.h)
        y = Math.max(Math.min(y, this.y_dim * this.h), this.h)

        let dx = 0.0
        let dy = 0.0

        let f

        switch (field) {
          case U_FIELD: f = this.u; dy = h2; break;
          case V_FIELD: f = this.v; dx = h2; break;
          case S_FIELD: f = this.d; dx = h2; dy = h2; break;
        }

        let x0 = Math.min(Math.floor((x - dx) * h1), this.x_dim - 1)
        let tx = ((x - dx) - x0 * this.h) * h1
        let x1 = Math.min(x0 + 1, this.x_dim - 1)

        let y0 = Math.min(Math.floor((y - dy) * h1), this.y_dim - 1)
        let ty = ((y - dy) - y0 * this.h) * h1
        let y1 = Math.min(y0 + 1, this.y_dim - 1)

        let sx = 1.0 - tx
        let sy = 1.0 - ty

        let val = sx * sy * f[x0 * n + y0] +
                  tx * sy * f[x1 * n + y0] +
                  tx * ty * f[x1 * n + y1] +
                  sx * ty * f[x0 * n + y1]

        return val
      }

      u_average (i, j) {
        const n = this.y_dim
        let u = (this.u[i*n + j-1] + this.u[i*n + j] +
                this.u[(i+1)*n + j-1] + this.u[(i+1)*n + j]) / 4
        
        return u
      }

      v_average (i, j) {
        const n = this.y_dim
        let v = (this.v[(i-1)*n + j] + this.v[i*n + j] +
                  this.v[(i-1)*n + j+1] + this.v[i*n + j+1]) / 4
      
        return v
      }

      advect_velocity (dt) {
        this.u_new.set(this.u)
        this.v_new.set(this.v)

        const n = this.y_dim
        const h2 = this.h / 2

        for (let i = 1; i < this.x_dim; i++) {
          for (let j = 1; j < this.y_dim; j++) {
            _count++

            if (this.s[i*n + j] != 0.0 && this.s[(i-1)*n + j] != 0.0 && j < this.y_dim - 1) {
              let x = i * this.h
              let y = j * this.h + h2
              let u = this.u[i*n + j]
              let v = this.v_average(i, j)
              v = this.sample_field(x, y, V_FIELD)
              x = x - u * dt
              y = y - v * dt
              u = this.sample_field(x, y, U_FIELD)
              this.u_new[i*n + j] = u
            }

            if (this.s[i*n + j] != 0.0 && this.s[i*n + j-1] != 0.0 && i < this.x_dim - 1) {
              let x = i * this.h + h2
              let y = j * this.h
              let u = this.u_average(i, j)
              u = this.sample_field(x, y, U_FIELD)
              let v = this.v[i*n + j]
              x = x - u * dt
              y = y - v * dt
              v = this.sample_field(x, y, V_FIELD)
              this.v_new[i*n + j] = v
            }
          }
        }
        this.u.set(this.u_new)
        this.v.set(this.v_new)
      }

      advect_smoke (dt) {
        this.d_new.set(this.d)

        const n = this.y_dim
        const h2 = this.h / 2

        for (let i = 1; i < this.x_dim - 1; i++) {
          for (let j = 1; j < this.y_dim - 1; j++) {
            if (this.s[i*n + j] != 0.0) {
              let u = (this.u[i*n + j] + this.u[(i+1)*n + j]) / 2
              let v = (this.v[i*n + j] + this.v[i*n + j+1]) / 2
              let x = i * this.h + h2 - u * dt
              let y = j * this.h + h2 - v * dt

              this.d_new[i*n + j] = this.sample_field(x, y, S_FIELD)
            }
          }
        }
        this.d.set(this.d_new)
      }

      simulate (dt, gravity, iterations) {
        this.integrate(dt, gravity)

        this.p.fill(0.0)
        this.solve_incompressibility(iterations, dt)

        this.extrapolate()
        this.advect_velocity(dt)
        this.advect_smoke(dt)
      }

      // fill cells with color gradient temp data
      test_fill () {
        let k = 0
        for (let i = 0; i < this.x_dim; i++) {
          for (let j = 0; j < this.y_dim; j++) {
            k++
            this.d[i * this.y_dim + j] = k / (this.x_dim * this.y_dim)
          }
        }
      }
    }


  return (
    <canvas id='fluid' className='absolute z-[-1] left-0 top-0 w-screen h-screen'></canvas>
  )
}

export default Sim