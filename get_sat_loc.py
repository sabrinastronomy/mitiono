#!/usr/bin/env python3
"""
Satellite Pass Query Script for Dominion Radio Astrophysical Observatory (DRAO)
Queries historical and future satellite passes with theta (elevation), phi (azimuth), and GPS data.

Location: DRAO, Penticton, BC, Canada
Coordinates: 49.3211°N, 119.6208°W, 545m elevation
"""

import requests
import json
from datetime import datetime, timedelta
import time
import math


class SatellitePassTracker:
    def __init__(self, api_key=None):
        """
        Initialize the satellite pass tracker.

        Args:
            api_key (str): N2YO API key (get from https://n2yo.com/api/)
        """
        self.api_key = api_key
        self.base_url = "https://api.n2yo.com/rest/v1/satellite"

        # DRAO coordinates
        self.drao_lat = 49.3211
        self.drao_lon = -119.6208
        self.drao_alt = 545  # meters

        # Common satellite NORAD IDs
        self.satellites = {
            'ISS': 25544,
            'HUBBLE': 20580,
            'AQUA': 27424,
            'TERRA': 25994,
            'LANDSAT_8': 39084,
            'LANDSAT_9': 49260,
            'SENTINEL_1A': 39634,
            'SENTINEL_1B': 41456,
            'SENTINEL_2A': 40697,
            'SENTINEL_2B': 42063,
            'NOAA_19': 33591,
            'NOAA_20': 43013
        }

    def get_satellite_passes(self, satellite_id, days=7, min_elevation=10):
        """
        Get visible satellite passes for DRAO location.

        Args:
            satellite_id (int): NORAD satellite ID
            days (int): Number of days to look ahead
            min_elevation (float): Minimum elevation angle in degrees

        Returns:
            dict: API response with pass data
        """
        if not self.api_key:
            print("Warning: No API key provided. Using demo data.")
            return self._get_demo_pass_data()

        url = f"{self.base_url}/visualpasses/{satellite_id}/{self.drao_lat}/{self.drao_lon}/{self.drao_alt}/{days}/{min_elevation}"

        params = {
            'apiKey': self.api_key
        }

        try:
            response = requests.get(url, params=params)
            response.raise_for_status()
            return response.json()
        except requests.exceptions.RequestException as e:
            print(f"Error fetching data: {e}")
            return None

    def get_satellite_positions(self, satellite_id, positions=10):
        """
        Get current satellite positions with azimuth and elevation.

        Args:
            satellite_id (int): NORAD satellite ID
            positions (int): Number of future positions to calculate

        Returns:
            dict: API response with position data
        """
        if not self.api_key:
            print("Warning: No API key provided. Using demo data.")
            return self._get_demo_position_data()

        url = f"{self.base_url}/positions/{satellite_id}/{self.drao_lat}/{self.drao_lon}/{self.drao_alt}/{positions}"

        params = {
            'apiKey': self.api_key
        }

        try:
            response = requests.get(url, params=params)
            response.raise_for_status()
            return response.json()
        except requests.exceptions.RequestException as e:
            print(f"Error fetching positions: {e}")
            return None

    def calculate_historical_pass(self, satellite_id, target_date, duration_hours=24):
        """
        Calculate historical satellite pass data for a specific date.
        Note: This requires TLE data from the target date for accurate results.

        Args:
            satellite_id (int): NORAD satellite ID
            target_date (str): Date in format 'YYYY-MM-DD'
            duration_hours (int): Duration to calculate passes for

        Returns:
            dict: Calculated pass data
        """
        print(f"Calculating historical pass data for {target_date}")
        print("Note: For accurate historical data, you need TLE data from the target date.")

        # This is a simplified example - real historical calculation would need:
        # 1. Historical TLE data from the target date
        # 2. Orbital mechanics calculations (SGP4/SDP4)
        # 3. Libraries like python-sgp4 or skyfield

        return {
            'satellite_id': satellite_id,
            'target_date': target_date,
            'observer_lat': self.drao_lat,
            'observer_lon': self.drao_lon,
            'observer_alt': self.drao_alt,
            'note': 'Historical calculation requires TLE data from target date'
        }

    def format_pass_data(self, pass_data):
        """
        Format pass data for display with theta, phi, and GPS coordinates.

        Args:
            pass_data (dict): Raw pass data from API

        Returns:
            str: Formatted pass information
        """
        if not pass_data or 'passes' not in pass_data:
            return "No pass data available"

        output = []
        output.append(f"Satellite Passes for DRAO (φ={self.drao_lat}°N, λ={self.drao_lon}°W)")
        output.append("=" * 70)

        for i, pass_info in enumerate(pass_data['passes'], 1):
            output.append(f"\nPass {i}:")
            output.append(
                f"  Date/Time: {datetime.fromtimestamp(pass_info['startUTC']).strftime('%Y-%m-%d %H:%M:%S UTC')}")
            output.append(f"  Duration: {pass_info['duration']} seconds")
            output.append(f"  Max Elevation (θ): {pass_info['maxEl']}°")
            output.append(f"  AOS Azimuth (φ): {pass_info['startAz']}°")
            output.append(f"  Max Azimuth (φ): {pass_info['maxAz']}°")
            output.append(f"  LOS Azimuth (φ): {pass_info['endAz']}°")
            output.append(f"  Visibility: {pass_info['mag']} mag")

        return "\n".join(output)

    def format_position_data(self, position_data):
        """
        Format position data showing theta, phi, and GPS coordinates.

        Args:
            position_data (dict): Raw position data from API

        Returns:
            str: Formatted position information
        """
        if not position_data or 'positions' not in position_data:
            return "No position data available"

        output = []
        output.append(f"Satellite Positions from DRAO (φ={self.drao_lat}°N, λ={self.drao_lon}°W)")
        output.append("=" * 70)

        for i, pos in enumerate(position_data['positions'], 1):
            timestamp = datetime.fromtimestamp(pos['timestamp'])
            output.append(f"\nPosition {i}:")
            output.append(f"  Time: {timestamp.strftime('%Y-%m-%d %H:%M:%S UTC')}")
            output.append(f"  Satellite GPS: φ={pos['satlatitude']}°, λ={pos['satlongitude']}°")
            output.append(f"  Altitude: {pos['sataltitude']} km")
            output.append(f"  Elevation (θ): {pos['elevation']}°")
            output.append(f"  Azimuth (φ): {pos['azimuth']}°")
            output.append(f"  Range: {pos['ra']} km")
            output.append(f"  Eclipsed: {'Yes' if pos['eclipsed'] else 'No'}")

        return "\n".join(output)

    def _get_demo_pass_data(self):
        """Return demo pass data when no API key is provided."""
        return {
            'passes': [
                {
                    'startUTC': int(time.time()) + 3600,
                    'duration': 420,
                    'maxEl': 45,
                    'startAz': 315,
                    'maxAz': 180,
                    'endAz': 45,
                    'mag': -2.5
                },
                {
                    'startUTC': int(time.time()) + 7200,
                    'duration': 380,
                    'maxEl': 30,
                    'startAz': 270,
                    'maxAz': 135,
                    'endAz': 90,
                    'mag': -1.8
                }
            ]
        }

    def _get_demo_position_data(self):
        """Return demo position data when no API key is provided."""
        current_time = int(time.time())
        return {
            'positions': [
                {
                    'timestamp': current_time,
                    'satlatitude': 51.5,
                    'satlongitude': -125.3,
                    'sataltitude': 408,
                    'elevation': 25,
                    'azimuth': 180,
                    'ra': 1200,
                    'eclipsed': False
                },
                {
                    'timestamp': current_time + 60,
                    'satlatitude': 52.1,
                    'satlongitude': -124.8,
                    'sataltitude': 410,
                    'elevation': 28,
                    'azimuth': 185,
                    'ra': 1150,
                    'eclipsed': False
                }
            ]
        }


def main():
    """Main function to demonstrate the satellite tracking functionality."""
    print("DRAO Satellite Pass Tracker")
    print("=" * 40)

    # Initialize tracker (replace 'YOUR_API_KEY' with actual N2YO API key)
    tracker = SatellitePassTracker(api_key=None)  # Set to None for demo

    print("Available satellites:")
    for name, norad_id in tracker.satellites.items():
        print(f"  {name}: {norad_id}")

    # Example: Get ISS passes
    print("\n1. Getting ISS passes for next 7 days:")
    iss_passes = tracker.get_satellite_passes(tracker.satellites['ISS'], days=7)
    print(tracker.format_pass_data(iss_passes))

    # Example: Get ISS current positions
    print("\n2. Getting ISS current positions:")
    iss_positions = tracker.get_satellite_positions(tracker.satellites['ISS'], positions=5)
    print(tracker.format_position_data(iss_positions))

    # Example: Historical calculation (requires additional implementation)
    print("\n3. Historical pass calculation:")
    historical = tracker.calculate_historical_pass(tracker.satellites['ISS'], '2024-01-15')
    print(json.dumps(historical, indent=2))

    print("\n" + "=" * 40)
    print("To use with real data:")
    print("1. Get API key from https://n2yo.com/api/")
    print("2. Replace 'YOUR_API_KEY' in the script")
    print("3. For historical data, consider using libraries like:")
    print("   - python-sgp4 for orbital calculations")
    print("   - skyfield for astronomical calculations")
    print("   - celestrak.com for historical TLE data")


if __name__ == "__main__":
    main()