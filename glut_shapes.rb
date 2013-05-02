require 'opengl'
require 'matrix'

#
# Methods translate radian/degree.
#
class Numeric
  def to_rad
    return self*Math::PI/180
  end
  def to_deg
    return self*180/Math::PI
  end
end

#
# Vector class extension to make easier drawing shapes with ruby-opengl.
#
class Vector
  include Math
  #
  # Method returns outer product of self and Vector o
  # ref) http://d.hatena.ne.jp/tanku/20100318/1268930228 (Japanese page)
  #
  def outer_product(o)
    raise 'vector size must be 3. (size#{self.size})' unless size == 3
    Vector[
      self[1] * o[2] - self[2] * o[1],
      self[2] * o[0] - self[0] * o[2],
      self[0] * o[1] - self[1] * o[0],
    ]
  end
  #
  # create unit vector[u,v,w] from self.
  # w = self/self.r (w and self have same direction.)
  #
  def create_unit_as_z
    w = self/self.r
    u0 = Vector[w[1]-w[2],  w[2]-w[0], w[0]-w[1]]
    u = u0/u0.r
    # Rotate u 90 degree around w.
    # ref: v = w(w・u) + (u-w(w・u))cos(90) - (u✕w)sin(90)
    v = w.inner_product(u)*w + u.outer_product(w)
    return [u, v, w]
  end
end


#
# GLUT class extension to make easier drawing shapes with ruby-opengl.
# Four shapes.
# -Ellipsoid
# -Macaroni
# -Cocoon
# -Legoblock

module GLUT
  extend Math
  # 
  # Draw Ellipsoid
  # 
  # ref) http://www.gamedev.net/topic/126624-generating-an-ellipsoid-in-opengl/
  def GLUT.Ellipsoid(x,y,z,slice,stack)
    tStep = Math::PI/slice.to_f
    sStep = Math::PI/stack.to_f
    (-Math::PI/2).step(Math::PI/2+0.0001, tStep){|t|
      GL.Begin(GL::TRIANGLE_STRIP)
      (-Math::PI).step(Math::PI+0.0001, sStep){|s|
        yield(t,s) if block_given?
        GL.Normal(cos(t)*cos(s), cos(t)*sin(s), sin(t))
        GL.Vertex(x*cos(t)*cos(s), y*cos(t)*sin(s), z*sin(t))
        GL.Vertex(x*cos(t+tStep)*cos(s), y*cos(t+tStep)*sin(s), z*sin(t+tStep))
      }
      GL.End
    }
  end
  #
  # Draw bended tube.
  # It is NOT guarantee tube length is len.
  # Give a block adequately.
  def GLUT.Macaroni(len, radius, slice, stack)
    return unless block_given?
    lStep = len/slice.to_f
    tStep = Math::PI/stack.to_f
    0.step(len-lStep*0.0001, lStep){|l|
      #p1,p2: Vectors discribe center of disks
      #u,v,w: Unit vector discribes coodinate in which draw disks.
      p1 = yield(l)
      p2 = yield(l+lStep)
      u1, v1, w1 = (p2-p1).create_unit_as_z
      trans_mat1 = Matrix.columns([u1,v1,w1])
      p3 = yield(l+2*lStep)
      u2, v2, w2 = (p3-p2).create_unit_as_z
      trans_mat2 = Matrix.columns([u2,v2,w2])
      GL.Begin(GL::TRIANGLE_STRIP)
      (-Math::PI).step(Math::PI+0.0001, tStep){|theta|
         GL.Normal(trans_mat1*Vector[cos(theta), sin(theta), 0])
         tr = trans_mat1 * Vector[radius*cos(theta), radius*sin(theta), 0]
         p0 = p1 + tr
         GL.Vertex(p0.to_a)
         p0 = p2 + trans_mat2 * Vector[radius*cos(theta+tStep), radius*sin(theta+tStep), 0]
         GL.Vertex(p0.to_a)
      }
      GL.End
    }
  end
  #
  # Draw cocoon (solid of revolution with changeable radius)
  # It is useful to draw human leg, cocoon, etc.
  # 
  def GLUT.Cocoon(height, slice, stack)
    return unless block_given?
    zStep = height/slice.to_f
    tStep = Math::PI/stack
    0.step(height-zStep+0.0001, zStep){|z|
       GL.Begin(GL::QUAD_STRIP)
      (-Math::PI).step(Math::PI+0.0001, tStep){|theta|
         r1 = yield(z, theta, tStep)
         r2 = yield(z+zStep, theta+tStep, tStep)
         phi = atan(zStep/(r2-r1).abs)
         GL.Normal(cos(theta)*sin(phi), sin(theta)*sin(phi), cos(phi))
         GL.Vertex(r1*cos(theta), r1*sin(theta), z)
         GL.Vertex(r2*cos(theta+tStep), r2*sin(theta+tStep), z+zStep)
      }
      GL.End
    }
  end
  #
  # Draw Lego BLock
  # Its size ratio is accurate.
  # 
def GLUT.Legoblock(size, size_unit=0.5)
    width_n = size[0]
    length_n = size[1]

    disks = Array.new(length_n*width_n)
    disks.map!{|e| e = GLU.NewQuadric() }
    tubes = Array.new(length_n*width_n + (length_n-1)*(width_n-1))
    tubes.map!{|e| e = GLU.NewQuadric() }

    w = size_unit*width_n  # Block width
    l = size_unit*length_n # Block length
    h = 1.2*size_unit      # Block height
    r = 0.5*5.0*size_unit/8.0 #Potch radius
    ph = 1.7*size_unit/8.0 #Potch height
    upr =sqrt(2)*size_unit - sqrt(2)*size_unit/2.0 - 0.5*5.0*size_unit/8.0
    
    GL.Begin(GL::QUADS)
      GL.Normal(0.0, 0.0, 1.0)
      GL.Vertex( w/2,  l/2, h/2)
      GL.Vertex( w/2, -l/2, h/2)
      GL.Vertex(-w/2, -l/2, h/2)
      GL.Vertex(-w/2,  l/2, h/2)
    GL.End
    
    GL.Begin(GL::QUADS)
      GL.Normal(1.0, 0.0, 0.0)
      GL.Vertex(w/2,  l/2,  h/2)
      GL.Vertex(w/2, -l/2,  h/2)
      GL.Vertex(w/2, -l/2, -h/2)
      GL.Vertex(w/2,  l/2, -h/2)
    GL.End
    GL.Begin(GL::QUADS)
      GL.Normal(-1.0, 0.0, 0.0)
      GL.Vertex(-w/2,  l/2,  h/2)
      GL.Vertex(-w/2, -l/2,  h/2)
      GL.Vertex(-w/2, -l/2, -h/2)
      GL.Vertex(-w/2,  l/2, -h/2)
    GL.End
    GL.Begin(GL::QUADS)
      GL.Normal(0.0, 1.0, 0.0)
      GL.Vertex( w/2, l/2, h/2)
      GL.Vertex(-w/2, l/2, h/2)
      GL.Vertex(-w/2, l/2, -h/2)
      GL.Vertex( w/2, l/2, -h/2)
    GL.End
    GL.Begin(GL::QUADS)
      GL.Normal(0.0, -1.0, 0.0)
      GL.Vertex( w/2, -l/2,  h/2)
      GL.Vertex(-w/2, -l/2,  h/2)
      GL.Vertex(-w/2, -l/2, -h/2)
      GL.Vertex( w/2, -l/2, -h/2)
    GL.End
    
    GL.PushMatrix()
    GL.Translate(0,0,h/2)
    i = 0
    dl = l.to_f/length_n
    dw = w.to_f/width_n
    offset_l = size_unit*l/length_n
    offset_w = size_unit*w/width_n
    width_n.times{|wi|
      length_n.times{|li|
        GL.PushMatrix()
        GL.Translate(offset_w+wi*dw-w.to_f/2, offset_l+li*dl-l.to_f/2, 0)
        GLU.Cylinder(tubes[i], r, r, ph, 20, 20)
        GL.Translate(0,0,ph)
        GL.Normal(0.0, 0.0, 1.0)
        GLU.Disk(disks[i], 0, r, 20, 20)
        GL.PopMatrix()
        i += 1
      }
    }
    GL.PopMatrix()
    
    dl = l.to_f/(length_n-1)
    dw = w.to_f/(width_n-1)
    offset_l = size_unit*l/(length_n-1)
    offset_w = size_unit*w/(width_n-1)

    GL.PushMatrix()
    GL.Translate(0,0,-h/2)
      (width_n-1).times{|wi|
      (length_n-1).times{|li|
        GL.PushMatrix()
        GL.Translate(offset_w+wi*dw-w.to_f/2, offset_l+li*dl-l.to_f/2, 0)
        GLU.Cylinder(tubes[i], upr, upr, h, 20, 20)
        GL.PopMatrix()
        i += 1
      }
    }
    GL.PopMatrix()
  end
end

