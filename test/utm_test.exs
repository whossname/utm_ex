defmodule UTMTest do
  use ExUnit.Case
  doctest UTM

  use PropCheck

  property "symetrical" do
    forall {lat, lon} <- {float(-80.0, 84.0), float(-180.0, 180.0)} do
      %{east: east, north: north} = UTM.from_wgs84(lat, lon)

      hemisphere = if lat >= 0, do: :north, else: :south
      zone = Integer.mod(floor((lon + 180) / 6), 60) + 1
      out = UTM.to_wgs84(east, north, zone, hemisphere)

      tolerance = 2 * :math.pow(10, -7)
      dlat = abs(lat - out.lat)
      dlon = abs(lon - out.lon)

      dlat < tolerance && dlon < tolerance
    end
  end
end
