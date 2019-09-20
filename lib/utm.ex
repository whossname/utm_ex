defmodule UTM do
  import :math

  @a 6_378_137.0
  @f 1 / 298.2572236
  @drad pi() / 180
  @k_0 0.9996
  @b @a * (1 - @f)
  @e sqrt(1 - pow(@b / @a, 2))
  @esq pow(@e, 2)
  @e0sq @esq / (1 - @esq)
  @e_1 (1 - sqrt(1 - @esq)) / (1 + sqrt(1 - @esq))

  @doc """
  Converts from WGS84 to UTM

  ## Examples:

      iex> UTM.from_wgs84(59.805241567229885, 11.40618711509996)
      %{east: 634980.0000762338, north: 6632172.011406599}

      iex> UTM.from_wgs84(63.50614385818526, 9.200909996819014)
      %{east: 5.1e5, north: 7042000.012790729}

      iex> UTM.from_wgs84(-31.953512, 115.857048)
      %{east: 391_984.4643429378, north:  6_464_146.921846279}
  """
  @mc_1 1 - @esq * (1 / 4 + @esq * (3 / 64 + 5 / 256 * @esq))
  @mc_2 @esq * (3 / 8 + @esq * (3 / 32 + 45 / 1024 * @esq))
  @mc_3 pow(@esq, 2) * (15 / 256 + @esq * 45 / 1204)
  @mc_4 pow(@esq, 3) * 35 / 3072

  def from_wgs84(lat, lon) do
    phi = lat * @drad
    zcm = 3 + 6 * (utmz(lon) - 1) - 180

    n = @a / sqrt(1 - pow(@e * sin(phi), 2))
    t = pow(tan(phi), 2)
    c = @e0sq * pow(cos(phi), 2)
    a = (lon - zcm) * @drad * cos(phi)

    x =
      @k_0 * n * a *
        (1 +
           pow(a, 2) *
             ((1 - t + c) / 6 + pow(a, 2) * (5 - 18 * t + pow(t, 2) + 72 * c - 58 * @e0sq) / 120)) +
        500_000

    m = phi * @mc_1
    m = m - sin(2 * phi) * @mc_2
    m = m + sin(4 * phi) * @mc_3
    m = m - sin(6 * phi) * @mc_4
    m = m * @a

    y =
      @k_0 *
        (m +
           n * tan(phi) *
             (pow(a, 2) *
                (1 / 2 +
                   pow(a, 2) *
                     ((5 - t + 9 * c + 4 * pow(c, 2)) / 24 +
                        pow(a, 2) * (61 - 58 * t + pow(t, 2) + 600 * c - 330 * @e0sq) / 720))))

    if y < 0 do
      %{east: x, north: y + 10_000_000}
    else
      %{east: x, north: y}
    end
  end

  @doc """
  Converts from UTM to WGS84

  ## Examples:

      iex> UTM.to_wgs84(634980.0, 6632172.0, 32, :north)
      %{lat: 59.805241567229885, lon: 11.40618711509996}

      iex> UTM.to_wgs84(510000, 7042000, 32, :north)
      %{lat: 63.50614385818526, lon: 9.200909996819014}

      iex> UTM.to_wgs84(391_984.4643429378, 6_464_146.921846279, 50, :south)
      %{lat: -31.95351191012423, lon: 115.8570480011093}
  """
  def to_wgs84(e, n, zone, hemisphere) do
    n =
      case hemisphere do
        :south -> n - 10_000_000
        :north -> n
      end

    m = n / @k_0

    mu = m / (@a * @mc_1)
    phi_1 = mu + @e_1 * (3 / 2 - 27 / 32 * pow(@e_1, 2)) * sin(2 * mu)
    phi_1 = phi_1 + pow(@e_1, 2) * (21 / 16 - 55 / 32 * pow(@e_1, 2)) * sin(4 * mu)
    phi_1 = phi_1 + pow(@e_1, 3) * (sin(6 * mu) * 151 / 96 + @e_1 * sin(8 * mu) * 1097 / 512)

    c_1 = @e0sq * pow(cos(phi_1), 2)
    t_1 = pow(tan(phi_1), 2)
    n_1 = @a / sqrt(1 - pow(@e * sin(phi_1), 2))
    r_1 = n_1 * (1 - @esq) / (1 - pow(@e * sin(phi_1), 2))
    d = (e - 500_000) / (n_1 * @k_0)

    phi =
      pow(d, 2) *
        (1 / 2 - pow(d, 2) * (5 + 3 * t_1 + 10 * c_1 - 4 * pow(c_1, 2) - 9 * @e0sq) / 24)

    phi =
      phi +
        pow(d, 6) * (61 + 90 * t_1 + 298 * c_1 + 45 * pow(t_1, 2) - 252 * @e0sq - 3 * pow(c_1, 2)) /
          720

    phi = phi_1 - n_1 * tan(phi_1) / r_1 * phi

    lon =
      d *
        (1 +
           pow(d, 2) *
             ((-1 - 2 * t_1 - c_1) / 6 +
                pow(d, 2) *
                  (5 - 2 * c_1 + 28 * t_1 - 3 * pow(c_1, 2) + 8 * @e0sq + 24 * pow(t_1, 2)) / 120)) /
        cos(phi_1)

    zcm = 3 + 6 * (zone - 1) - 180
    %{lat: phi / @drad, lon: zcm + lon / @drad}
  end

  defp utmz(lon), do: 1 + Float.floor((lon + 180) / 6)
end
