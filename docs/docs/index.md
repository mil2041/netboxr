### GSoC22 Reports
<table>
  {% for post in site.reports reversed %}
    <tr>
      <td>{{ post.date | date_to_long_string }}</td>
      <td><a href="{{ post.url | remove_first:'/' }}">{{ post.title }}</a></td>
    </tr>
  {% endfor %}
</table>


### Contributors
<ul>
  {% for member in site.data.contributors %}
      <li>
        <a target="_blank" href="https://github.com/{{member.github}}">{{ member.name }} ({{ member.position }})</a>
      </li>
  {% endfor %}
</ul>
